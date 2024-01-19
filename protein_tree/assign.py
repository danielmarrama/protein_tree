#!/usr/bin/env python3

import warnings
warnings.filterwarnings('ignore')

import os
import glob
import pandas as pd

from ARC.classifier import SeqClassifier
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pathlib import Path
from pepmatch import Preprocessor, Matcher

from protein_tree.get_data import DataFetcher

class GeneAndProteinAssigner:
  def __init__(
    self,
    taxon_id,
    species_path,
    is_vertebrate,
    num_threads,
    build_path = Path(__file__).parent.parent / 'build',
    bin_path = Path(__file__).parent.parent / 'bin',
  ):
    
    self.species_path = species_path
    self.taxon_id = taxon_id
    self.is_vertebrate = is_vertebrate
    self.build_path = build_path
    self.bin_path = bin_path
    self.num_threads = num_threads

    # initialize dicts for assignments
    self.source_gene_assignment = {}
    self.source_protein_assignment = {}
    self.source_assignment_score = {}
    self.source_arc_assignment = {}
    self.epitope_protein_assignment = {}

    # create UniProt ID -> gene map
    self.proteome = pd.read_csv(f'{self.species_path}/proteome.tsv', sep='\t')
    self.uniprot_id_to_gene_map = dict(
      zip(
        self.proteome['Protein ID'], 
        self.proteome['Gene']
      )
    )
    # create UniProt ID -> protein name map
    self.uniprot_id_to_name_map = dict(
      zip(
        self.proteome['Protein ID'], 
        self.proteome['Protein Name']
      )
    )

  def assign(self, sources_df: pd.DataFrame, epitopes_df: pd.DataFrame) -> None:
    """Overall function to assign genes and parent proteins to sources and epitopes.

    Args:
      sources_df: DataFrame of source antigens for a species.
      epitopes_df: DataFrame of epitopes for a species.
    """
    # assign None to all source antigens and epitopes to start
    for i, row in sources_df.iterrows():
      self.source_gene_assignment[row['Accession']] = None
      self.source_protein_assignment[row['Accession']] = None
      self.source_assignment_score[row['Accession']] = None
    for i, row in epitopes_df.iterrows():
      self.epitope_protein_assignment[(row['Accession'], row['Sequence'])] = None

    self.source_to_epitopes_map = self._create_source_to_epitopes_map(epitopes_df)

    sources_df['Length'] = sources_df['Sequence'].str.len() # add length column to sources
    self.source_length_map = dict( # create map of source antigens to their length
      zip(
        sources_df['Accession'],
        sources_df['Length']
      )
    )
    print('Assigning source antigens...')
    self._assign_sources(sources_df)
    print('Done assigning source antigens.\n')

    self._assign_allergens()
    self._assign_manuals()

    print('Assigning epitopes...')
    self._assign_epitopes(epitopes_df)
    print('Done.\n')

    # map source antigens to their BLAST matches and ARC assignments
    sources_df.loc[:, 'Assigned Gene'] = sources_df['Accession'].map(self.source_gene_assignment)
    sources_df.loc[:, 'Assigned Protein ID'] = sources_df['Accession'].map(self.source_protein_assignment)
    sources_df.loc[:, 'Assigned Protein Name'] = sources_df['Assigned Protein ID'].map(self.uniprot_id_to_name_map)
    sources_df.loc[:, 'Assignment Score'] = sources_df['Accession'].map(self.source_assignment_score)
    sources_df.loc[:, 'ARC Assignment'] = sources_df['Accession'].map(self.source_arc_assignment)

    # map epitope source antigens to assignments above and then PEPMatch assignments
    epitopes_df.loc[:, 'Assigned Gene'] = epitopes_df['Accession'].map(self.source_gene_assignment)
    epitopes_df.set_index(['Accession', 'Sequence'], inplace=True)
    epitopes_df.loc[:, 'Assigned Protein ID'] = epitopes_df.index.map(self.epitope_protein_assignment)
    epitopes_df.loc[:, 'Assigned Protein Name'] = epitopes_df['Assigned Protein ID'].map(self.uniprot_id_to_name_map)
    epitopes_df.reset_index(inplace=True)
    epitopes_df.loc[:, 'ARC Assignment'] = epitopes_df['Accession'].map(self.source_arc_assignment)

    epitopes_df.drop_duplicates(subset=['Accession', 'Sequence'], inplace=True) # drop duplicate epitopes
    sources_df.drop(columns=['Sequence'], inplace=True) # drop sequence column for output

    self._remove_files()
    
    num_sources = len(sources_df['Accession'].drop_duplicates())
    num_epitopes = len(epitopes_df[['Accession', 'Sequence']].drop_duplicates())
    num_matched_sources = len(sources_df[sources_df['Assigned Protein ID'].notnull()])
    num_matched_epitopes = len(epitopes_df[epitopes_df['Assigned Protein ID'].notnull()])
  
    assigner_data = (
      num_sources,
      num_epitopes,
      num_matched_sources,
      num_matched_epitopes
    )

    return assigner_data, epitopes_df, sources_df

  def _assign_sources(self, sources_df: pd.DataFrame) -> None:
    """Assign a gene to the source antigens of a species.

    Run ARC for vertebrates to assign MHC/TCR/Ig to source antigens first.
    Run BLAST for all other source antigens to assign a gene and protein. 
    If there are ties, use PEPMatch to search the epitopes within the protein 
    sequences andselect the protein with the most matches.

    Args:
      sources_df: DataFrame of source antigens for a species.
    """    
    self._write_to_fasta(sources_df, 'sources') # write sources to FASTA file

    print('Running BLAST for source antigens...')
    self._run_blast()

    if self.is_vertebrate:

      if sources_df['Sequence'].isna().all():
        return
      
      print('Running ARC for MHC/TCR/Ig assignments...')
      self._run_arc(sources_df)

  def _create_source_to_epitopes_map(self, epitopes_df: pd.DataFrame) -> dict:
    """Create a map from source antigens to their epitopes.
    
    Args:
      epitopes_df: DataFrame of epitopes for a species.
    """    
    source_to_epitopes_map = {}
    for i, row in epitopes_df.iterrows():
      if row['Accession'] in source_to_epitopes_map.keys():
        source_to_epitopes_map[row['Accession']].append(row['Sequence'])
      else:
        source_to_epitopes_map[row['Accession']] = [row['Sequence']]
    
    return source_to_epitopes_map 

  def _write_to_fasta(self, sources_df: pd.DataFrame, filename: str) -> None:
    """Write source antigens to FASTA file.

    Args:
      sources_df: DataFrame of source antigens for a species.
      filename: Name of the FASTA file to write to.    
    """          
    seq_records = [] # create seq records of sources with ID and sequence
    for i, row in sources_df.iterrows():

      # if there is no sequence, use empty string
      if pd.isnull(row['Sequence']):
        continue

      seq_records.append(
        SeqRecord(
          Seq(row['Sequence']),
          id=row['Accession'],
          description='')
      )
    with open(f'{self.species_path}/{filename}.fasta', 'w') as f:
      SeqIO.write(seq_records, f, 'fasta')

  def _preprocess_proteome_if_needed(self) -> None:
    """Preprocess the proteome if the preprocessed files don't exist."""
    if not os.path.exists(f'{self.species_path}/proteome.db'):
      gp_proteome = f'{self.species_path}/gp_proteome.fasta' if os.path.exists(f'{self.species_path}/gp_proteome.fasta') else ''
      Preprocessor(
        proteome = f'{self.species_path}/proteome.fasta',
        preprocessed_files_path = f'{self.species_path}',
        gene_priority_proteome=gp_proteome
      ).sql_proteome(k = 5)

  def _run_blast(self) -> None:
    """BLAST source antigens against the selected proteome, then read in with
    pandas and assign column names. By default, blastp doesn't return header.

    Then, create a quality score based on % identity, alignment length, and
    query length. Select the best match for each source antigen and assign
    the gene symbol and UniProt ID to the source_gene_assignment and
    source_protein_assignment maps.

    If there are ties, use PEPMatch to search the epitopes within the protein
    sequences and select the protein with the most matches.
    """
    # escape parentheses in species path
    species_path = str(self.species_path).replace('(', '\\(').replace(')', '\\)')

    os.system( # make BLAST database from proteome
      f'{self.bin_path}/makeblastdb -in {species_path}/proteome.fasta '\
      f'-dbtype prot > /dev/null'
    )
    os.system( # run blastp
      f'{self.bin_path}/blastp -query {species_path}/sources.fasta '\
      f'-db {species_path}/proteome.fasta '\
      f'-evalue 1  -num_threads {self.num_threads} -outfmt 10 '\
      f'-out {species_path}/blast_results.csv'
    )
    result_columns = [
      'Query', 'Target', 'Percentage Identity', 'Alignment Length', 
      'Mismatches', 'Gap Opens', 'Query Start', 'Query End', 
      'Target Start', 'Target End', 'e-Value', 'Bit Score'
    ]
    blast_results_df = pd.read_csv( # read in BLAST results
      f'{self.species_path}/blast_results.csv', names=result_columns
    )

    # extract the UniProt ID from the target column
    blast_results_df['Target'] = blast_results_df['Target'].str.split('|').str[1]

    # take out "-#" portion of the target UniProt ID because these are isoforms and 
    # won't be mapped properly to gene symbols
    blast_results_df['Target'] = blast_results_df['Target'].str.split('-').str[0]
    
    # map target UniProt IDs to gene symbols
    blast_results_df['Target Gene Symbol'] = blast_results_df['Target'].map(self.uniprot_id_to_gene_map)

    # create a quality score based on % identity, alignment length, and query length
    blast_results_df['Query Length'] = blast_results_df['Query'].astype(str).map(self.source_length_map)
    blast_results_df['Quality Score'] = blast_results_df['Percentage Identity'] * (blast_results_df['Alignment Length'] / blast_results_df['Query Length'])

    # join proteome metadata to select the best match for each source antigen
    blast_results_df = blast_results_df.merge(
      self.proteome[['Protein ID', 'Protein Existence Level', 'Gene Priority']], 
      left_on='Target', 
      right_on='Protein ID', 
      how='left'
    )
    # sort by quality score (descending), gene priority (descending), and
    # protein existence level (ascending)
    blast_results_df.sort_values(
      by=['Quality Score', 'Gene Priority', 'Protein Existence Level'],
      ascending=[False, False, True],
      inplace=True
    )
    # after sorting, drop duplicates based on 'Query', keeping only the first (i.e., best) match.
    blast_results_df.drop_duplicates(subset='Query', keep='first', inplace=True)

    # assign gene symbols, protein ID, and score
    for i, row in blast_results_df.iterrows():
      self.source_gene_assignment[str(row['Query'])] = row['Target Gene Symbol']
      self.source_protein_assignment[str(row['Query'])] = row['Target']
      self.source_assignment_score[str(row['Query'])] = row['Quality Score']

  def _assign_epitopes(self, epitopes_df: pd.DataFrame) -> None:
    """Assign a parent protein to each epitope.
    
    Preprocess the proteome and then search all the epitopes within
    the proteome using PEPMatch. Then, assign the parent protein
    to each epitope by selecting the best isoform of the assigned gene for
    its source antigen.
    """
    self._preprocess_proteome_if_needed()

    # search all epitopes within the proteome using PEPMatch
    all_epitopes = epitopes_df['Sequence'].unique().tolist()
    all_matches_df = Matcher(
      query = all_epitopes,
      proteome_file = f'{self.species_path}/proteome.fasta',
      max_mismatches = 0, 
      k = 5, 
      preprocessed_files_path = f'{self.species_path}', 
      best_match=False, 
      output_format='dataframe',
      sequence_version=False
    ).match()
    
    # if no source antigens were assigned, return
    if not self.source_gene_assignment or not self.source_protein_assignment:
      return

    # create dataframes of source antigen mappings so we can merge and perform operations
    epitope_source_map_df = pd.DataFrame({
      'Epitope': epitope,
      'Source Antigen': antigen
    } for antigen, epitopes in self.source_to_epitopes_map.items() for epitope in epitopes)
    source_protein_assignment_df = pd.DataFrame({
      'Source Antigen': antigen, 
      'Protein ID': protein_id
    } for antigen, protein_id in self.source_protein_assignment.items())
    source_gene_assignment_df = pd.DataFrame({
      'Source Antigen': antigen, 
      'Gene': gene
    } for antigen, gene in self.source_gene_assignment.items())

    # merge the source antigen for each epitope
    merged_df = pd.merge( 
      all_matches_df, 
      epitope_source_map_df, 
      how='left', 
      left_on='Query Sequence', 
      right_on='Epitope'
    )
    # merge the protein ID for each source antigen
    merged_df = pd.merge( 
      merged_df, 
      source_protein_assignment_df, 
      how='left', 
      on='Source Antigen',
      suffixes=('', '_assigned')
    )

    # now, grab epitopes that match to the assigned protein for their source antigen
    assigned_epitopes = merged_df[merged_df['Protein ID'] == merged_df['Protein ID_assigned']]

    self.epitope_protein_assignment.update(
      dict(zip(
        zip(assigned_epitopes['Source Antigen'], assigned_epitopes['Query Sequence']),
        assigned_epitopes['Protein ID']))
    )
    # remove assigned epitopes
    merged_df = merged_df[~merged_df['Query Sequence'].isin(assigned_epitopes['Query Sequence'])]

    # grab the rest of the epitopes and merged with gene assignments
    unassigned_epitopes = merged_df[merged_df['Protein ID'] != merged_df['Protein ID_assigned']]

    if unassigned_epitopes.empty: # if there are no unassigned epitopes, return
      return len(self.epitope_protein_assignment.keys())

    # merge with gene assignments
    unassigned_epitopes = pd.merge(
      unassigned_epitopes,
      source_gene_assignment_df,
      how='left',
      on='Source Antigen',
      suffixes=('', '_assigned')
    )
    # now, get the isoform of the assigned gene with the best protein existence level
    best_isoform_indices = unassigned_epitopes.groupby(['Gene', 'Query Sequence'])['Protein Existence Level'].idxmin()
    best_isoforms = unassigned_epitopes.loc[best_isoform_indices, ['Gene', 'Query Sequence', 'Protein ID']].set_index(['Gene', 'Query Sequence'])['Protein ID']
    unassigned_epitopes['Best Isoform ID'] = unassigned_epitopes.set_index(['Gene', 'Query Sequence']).index.map(best_isoforms)
    
    # drop any unassigned epitopes that couldn't be assigned an isoform
    unassigned_epitopes = unassigned_epitopes.dropna(subset=['Best Isoform ID'])

    # Update the protein assignment with the best isoform
    self.epitope_protein_assignment.update(
      dict(zip(
        zip(unassigned_epitopes['Source Antigen'], unassigned_epitopes['Query Sequence']),
        unassigned_epitopes['Best Isoform ID']))
    )

  def _run_arc(self, sources_df: pd.DataFrame) -> None:
    """Run ARC to assign MHC/TCR/Ig to source antigens.
    
    Args:
      sources_df: DataFrame of source antigens for a species.
    """
    try: # read ARC results if they already exist and only run ARC on new sources
      past_arc_results_df = pd.read_csv(f'{self.species_path}/ARC_results.tsv', sep='\t')
      past_arc_ids = past_arc_results_df['id'].tolist()
      arc_sources_df = sources_df[~sources_df['Accession'].isin(past_arc_ids)]
      self._write_to_fasta(arc_sources_df, 'ARC_sources')
      past_results = True
    except FileNotFoundError:
      self._write_to_fasta(sources_df, 'ARC_sources')
      past_results = False

    # make sure ARC_sources.fasta isn't empty, then run ARC
    if os.stat(f'{self.species_path}/ARC_sources.fasta').st_size != 0:
      SeqClassifier( # run ARC
        outfile = f'{self.species_path}/ARC_results.tsv',
        threads = self.num_threads,
        hmmer_path = str(self.bin_path) + '/',
        blast_path = str(self.bin_path) + '/',
      ).classify_seqfile(f'{self.species_path}/ARC_sources.fasta')

    os.remove(f'{self.species_path}/ARC_sources.fasta')
    
    arc_results_df = pd.read_csv(f'{self.species_path}/ARC_results.tsv', sep='\t')
    if past_results: # combine old results and new, then write all to file
      arc_results_df = pd.concat([arc_results_df, past_arc_results_df])
      arc_results_df.to_csv(f'{self.species_path}/ARC_results.tsv', sep='\t', index=False)

    if not arc_results_df.dropna(subset=['class']).empty:
      self.source_arc_assignment = arc_results_df.set_index('id')['class'].to_dict()

  def _assign_allergens(self) -> None:
    """Get allergen data from allergen.org and then assign allergens to sources."""
    allergen_df = pd.read_csv(self.build_path / 'arborist' / 'allergens.tsv', sep='\t')
    allergen_map = allergen_df.set_index('AccProtein')['Name'].to_dict()

    for k, v in self.source_protein_assignment.items():
      if v in allergen_map.keys():
        self.uniprot_id_to_name_map[v] = allergen_map[v]

  def _assign_manuals(self) -> None:
    """Get manual assignments from manual_assignments.tsv and then assign
    genes and proteins to sources.
    """
    # manual_assignments.tsv should be in the directory above this one
    manual_df = pd.read_csv(self.build_path / 'arborist' / 'manual_assignments.tsv', sep='\t')

    manual_gene_map = manual_df.set_index('Accession')['Accession Gene'].to_dict()
    manual_protein_id_map = manual_df.set_index('Accession')['Parent Accession'].to_dict()
    manual_protein_name_map = manual_df.set_index('Accession')['Parent Name'].to_dict()
    
    for k, v in self.source_gene_assignment.items():
      if k in manual_gene_map.keys():
        self.source_gene_assignment[k] = manual_gene_map[k]
        self.source_protein_assignment[k] = manual_protein_id_map[k]
        self.source_assignment_score[k] = -1
        self.uniprot_id_to_name_map[k] = manual_protein_name_map[k]

  def _remove_files(self) -> None:
    """Delete all the files that were created when making the BLAST database."""
    for extension in ['pdb', 'phr', 'pin', 'pjs', 'pot', 'psq', 'ptf', 'pto']:
      try: # if DB wasn't created this will throw an error
        os.remove(glob.glob(f'{self.species_path}/*.{extension}')[0])
      except IndexError:
        pass 
    
    os.remove(f'{self.species_path}/blast_results.csv')
    os.remove(f'{self.species_path}/sources.fasta')

    # try: # if proteome.db exists, remove it
    #   os.remove(f'{self.species_path}/proteome.db')
    # except OSError:
    #   pass

def run(taxon_id, species_name, group, all_taxa, build_path, bin_path, all_epitopes, all_antigens, num_threads):
  species_path = build_path / 'species' / f'{taxon_id}' # directory to write species data

  epitopes_df = DataFetcher(build_path).get_epitopes_for_species(all_epitopes, all_taxa)
  sources_df = DataFetcher(build_path).get_sources_for_species(all_antigens, epitopes_df['Accession'].tolist())

  if sources_df.empty or epitopes_df.empty:
    return

  is_vertebrate = group == 'vertebrate'

  print(f'Assigning peptides and sources for species w/ taxon ID: {taxon_id}).')
  Assigner = GeneAndProteinAssigner(
    taxon_id,
    species_path,
    is_vertebrate,
    num_threads=num_threads,
    build_path=build_path,
    bin_path=bin_path
  )

  assigner_data, epitope_assignments, antigen_assignments = Assigner.assign(sources_df, epitopes_df)
  num_sources, num_epitopes, num_matched_sources, num_matched_epitopes = assigner_data

  epitope_assignments.to_csv(species_path / 'epitope-assignments.tsv', sep='\t', index=False)
  antigen_assignments.to_csv(species_path / 'antigen-assignments.tsv', sep='\t', index=False)

  assignment_data = (
    taxon_id,
    species_name,
    num_sources,
    num_epitopes,
    num_matched_sources,
    num_matched_epitopes
  )

  species_data = pd.read_csv(species_path / 'species-data.tsv', sep='\t')
  species_data.loc[:, 'Num Sources'] = num_sources
  species_data.loc[:, 'Num Peptides'] = num_epitopes
  species_data.loc[:, '% Assigned Sources'] = round((num_matched_sources / num_sources)*100, 2)
  species_data.loc[:, '% Assigned Peptides'] = round((num_matched_epitopes / num_epitopes)*100, 2)
  species_data.to_csv(species_path / 'species-data.tsv', sep='\t', index=False)

  return assignment_data

if __name__ == '__main__':
  import argparse

  parser = argparse.ArgumentParser()
  
  parser.add_argument(
    '-b', '--build_path', 
    type=str,
    default=Path(__file__).parent.parent / 'build',
    help='Path to data directory.'
  )
  parser.add_argument(
    '-B', '--bin_path',
    type=str,
    default=Path(__file__).parent.parent / 'bin',
    help='Path to bin directory.'
  )
  parser.add_argument(
    '-a', '--all_species', 
    action='store_true', 
    help='Build protein tree for all IEDB species.'
  )
  parser.add_argument(
    '-t', '--taxon_id', 
    type=int,
    help='Taxon ID for the species to pull data for.'
  )
  parser.add_argument(
    '-n', '--num_threads', 
    type=int,
    default=os.cpu_count() - 1,
    help='Number of threads to use.'
  )

  args = parser.parse_args()

  build_path = Path(args.build_path)
  bin_path = Path(args.bin_path)
  all_species = args.all_species
  taxon_id = args.taxon_id
  num_threads = args.num_threads

  assert all_species or taxon_id, 'Please specify either --all_species or --taxon_id.'
  assert (all_species and not taxon_id) or (taxon_id and not all_species), 'Please specify either --all_species or --taxon_id.'

  # TODO: replace the data/active-species.tsv with updated arborist active-species somehow
  species_df = pd.read_csv(build_path / 'arborist' / 'active-species.tsv', sep='\t')
  valid_taxon_ids = species_df['Species ID'].tolist()

  all_epitopes = DataFetcher(build_path).get_all_epitopes()
  all_antigens = DataFetcher(build_path).get_all_antigens()

  all_taxa_map = dict(zip( # map taxon ID to list of all children taxa
    species_df['Species ID'],
    species_df['Active Taxa']))
  taxon_to_species_map = dict(zip( # map taxon ID to species name
    species_df['Species ID'],
    species_df['Species Label']))

  if all_species: # run all species at once
    for taxon_id in valid_taxon_ids:
      species_name = taxon_to_species_map[taxon_id]
      group = species_df[species_df['Species ID'] == taxon_id]['Group'].iloc[0]
      all_taxa = [int(taxon) for taxon in all_taxa_map[taxon_id].split(', ')]
      assignment_data = run(
        taxon_id, species_name, group, all_taxa, build_path, bin_path, all_epitopes, all_antigens, num_threads
      )

  else: # one species at a time
    assert taxon_id in valid_taxon_ids, f'{taxon_id} is not a valid taxon ID.'
    species_name = taxon_to_species_map[taxon_id]
    all_taxa = [int(taxon) for taxon in all_taxa_map[taxon_id].split(', ')]
    group = species_df[species_df['Species ID'] == taxon_id]['Group'].iloc[0]
    assignment_data = run(
      taxon_id, species_name, group, all_taxa, build_path, bin_path, all_epitopes, all_antigens, num_threads
    )