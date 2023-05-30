#!/usr/bin/env python3

import os
import glob
import pandas as pd

from ARC.classifier import SeqClassifier
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pathlib import Path
from pepmatch import Preprocessor, Matcher


class GeneAndProteinAssigner:
  def __init__(self, taxon_id, species_name, is_vertebrate):
    self.taxon_id = taxon_id
    self.species_path = Path(f'species/{taxon_id}-{species_name.replace(" ", "_")}')
    self.is_vertebrate = is_vertebrate

    # create UniProt ID to gene symbol map from proteome.csv file
    proteome = pd.read_csv(f'{self.species_path}/proteome.csv')
    self.uniprot_id_to_gene_symbol_map = dict(
      zip(
        proteome['UniProt ID'], 
        proteome['Gene Symbol']
    ))
    self.source_gene_assignment = {}
    self.source_protein_assignment = {}


  def assign(self, sources_df: pd.DataFrame, epitopes_df: pd.DataFrame) -> None:
    """Overall function to assign genes and parent proteins to sources and epitopes.

    Args:
      sources_df: DataFrame of source antigens for a species.
      epitopes_df: DataFrame of epitopes for a species.
    """
    num_sources = len(sources_df['Accession'].unique())
    num_epitopes = len(epitopes_df['Sequence'].unique())
    
    # create source to epitope map
    self.source_to_epitopes_map = self._create_source_to_epitopes_map(epitopes_df)
    self.source_length_map = dict( # create map of source antigens to their length
      zip(
        sources_df['Accession'],
        sources_df['Length']
      )
    )

    num_matched_sources = self._assign_genes(sources_df, epitopes_df, num_sources)
    num_matched_epitopes = self._assign_parents(epitopes_df)

    self._assign_allergens()
    self._assign_manuals()

    # map source antigens to their best blast matches (UniProt ID and gene) for sources
    sources_df['Assigned Gene'] = sources_df['Accession'].map(self.source_gene_assignment)
    sources_df['Assigned Protein ID'] = sources_df['Accession'].map(self.source_protein_assignment)
    
    # map source antigens to their best blast matches (gene) for epitopes
    epitopes_df['Assigned Gene'] = epitopes_df['Source Accession'].map(self.source_gene_assignment)
    epitopes_df['Assigned Parent Protein ID'] = epitopes_df['Sequence'].map(self.epitope_parent_assignment)

    sources_df.drop(columns=['Sequence'], inplace=True) # drop sequence column for output
    sources_df.to_csv(f'{self.species_path}/sources.csv', index=False)
    epitopes_df.to_csv(f'{self.species_path}/epitopes.csv', index=False)

    self._remove_files()
    
    assigner_data = (
      num_sources,
      num_epitopes,
      num_matched_sources,
      num_matched_epitopes
    )
    return assigner_data


  def _assign_genes(
    self, sources_df: pd.DataFrame, epitopes_df: pd.DataFrame, num_sources: int
  ) -> None:
    """Assign a gene to the source antigens of a species.

    Run ARC for vertebrates to assign MHC/TCR/Ig to source antigens first.

    Then, run MMseqs2 if there are more than 1000 source antigens, for speed.
    Run BLAST for the remaining ones that don't have a match below e-value 10^-100.

    Otherwise, run BLAST to assign a gene to each source antigen. If there are
    ties, use PEPMatch to search the epitopes within the protein sequences and
    select the protein with the most matches.

    Args:
      sources_df: DataFrame of source antigens for a species.
      epitopes_df: DataFrame of epitopes for a species.
    """    
    self._sources_to_fasta(sources_df) # write sources to FASTA file

    if self.is_vertebrate:
      self._run_arc()

    if num_sources > 1000:
      self._run_mmseqs2()
      if self.source_gene_assignment.keys() != sources_df['Accession'].tolist():
        self._run_blast()
    else:
      self._run_blast()

    num_matched_sources = len(self.source_gene_assignment.keys())

    return num_matched_sources   


  def _create_source_to_epitopes_map(self, epitopes_df: pd.DataFrame) -> dict:
    """Create a map from source antigens to their epitopes.
    Args:
      epitopes_df: DataFrame of epitopes for a species.
    """    
    source_to_epitopes_map = {}
    for i, row in epitopes_df.iterrows():
      if row['Source Accession'] in source_to_epitopes_map.keys():
        source_to_epitopes_map[row['Source Accession']].append(row['Sequence'])
      else:
        source_to_epitopes_map[row['Source Accession']] = [row['Sequence']]
    
    return source_to_epitopes_map 


  def _sources_to_fasta(self, sources_df: pd.DataFrame) -> None:
    """Write source antigens to FASTA file."""          
    seq_records = [] # create seq records of sources with ID and sequence
    for i, row in sources_df.iterrows():
      seq_records.append(
        SeqRecord(
          Seq(row['Sequence']),
          id=row['Accession'],
          description='')
      )
    with open(f'{self.species_path}/sources.fasta', 'w') as f:
      SeqIO.write(seq_records, f, 'fasta')


  def _run_arc(self) -> None:
    """Run ARC to assign MHC/TCR/Ig to source antigens."""
    # pass sources.fasta to ARC
    SeqClassifier(
      outfile=f'{self.species_path}/ARC_results.tsv',
      blast_path = './' 
    ).classify_seqfile(f'{self.species_path}/sources.fasta')
    arc_results = pd.read_csv(f'{self.species_path}/ARC_results.tsv', sep='\t')
    
    if not arc_results.dropna(subset=['class']).empty:
      arc_assignments = arc_results.set_index('id')['class'].to_dict()
      self.source_gene_assignment.update(arc_assignments)


  def _run_mmseqs2(self) -> None:
    pass


  def _preprocess_proteome_if_needed(self) -> None:
    """Preprocess the proteome if the preprocessed files don't exist."""
    if not os.path.exists(f'{self.species_path}/proteome.db'):
        gp_proteome = f'{self.species_path}/gp_proteome.fasta' if os.path.exists(f'{self.species_path}/gp_proteome.fasta') else ''
        Preprocessor(
          proteome = f'{self.species_path}/proteome.fasta',
          preprocessed_files_path = f'{self.species_path}',
          gene_priority_proteome=gp_proteome
        ).sql_proteome(k = 5)


  def _search_epitopes(self, epitopes: list, best_match: bool = True) -> pd.DataFrame:
    """Search epitopes within the proteome using PEPMatch."""
    df = Matcher(
      query = epitopes,
      proteome_file = f'{self.species_path}/proteome.fasta',
      max_mismatches = 0, 
      k = 5, 
      preprocessed_files_path = f'{self.species_path}', 
      best_match=best_match, 
      output_format='dataframe'
    ).match()
    return df


  def _assign_parents(self, epitopes_df: pd.DataFrame) -> tuple:
    """Assign a parent protein to each epitope.
    
    Preprocess the proteome and then search all the epitopes within
    the proteome using PEPMatch. Then, assign the parent protein
    to each epitope by selecting the best isoform of the assigned gene for
    its source antigen.

    Args:
      epitopes_df: DataFrame of epitopes for a species.
    """
    self._preprocess_proteome_if_needed()

    # loop through the source antigens so we can limit the epitope matches
    # to the assigned gene for each source antigen
    all_matches_df = pd.DataFrame()
    for antigen, epitopes in self.source_to_epitopes_map.items():
      matches_df = self._search_epitopes(epitopes, best_match=True)

      # try to isolate matches to the assigned gene for the source antigen
      if antigen in self.source_gene_assignment.keys():
        matches_df = matches_df[matches_df['Gene'] == self.source_gene_assignment[antigen]]

      # if there aren't best matches that match the assigned gene, search again without best match
      if matches_df.empty:
        matches_df = self._search_epitopes(epitopes, best_match=False)

        # if there are ties, select the protein with the best protein existence level
        index = matches_df.groupby(['Query Sequence'])['Protein Existence Level'].transform(min) == matches_df['Protein Existence Level']
        matches_df = matches_df[index]
      
      # concatenate the matches to the all_matches_df
      all_matches_df = pd.concat([all_matches_df, matches_df])

    self.epitope_parent_assignment = dict(
      zip(all_matches_df['Query Sequence'], 
      all_matches_df['Protein ID'])
    )

    # count the number of epitopes with matches
    num_matched_epitopes = len(
      all_matches_df.dropna(subset=['Matched Sequence'])['Query Sequence'].unique()
    )

    return num_matched_epitopes


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
    species_path = str(self.species_path).replace('(', '\(').replace(')', '\)')

    os.system( # make BLAST database from proteome
      f'./makeblastdb -in {species_path}/proteome.fasta -dbtype prot'
    ) 
    os.system( # run blastp
      f'./blastp -query {species_path}/sources.fasta '\
      f'-db {species_path}/proteome.fasta '\
      f'-evalue 1 -num_threads 12 -outfmt 10 '\
      f'-out {species_path}/blast_results.csv'
    ) 
    result_columns = [
      'Query', 'Subject', 'Percentage Identity', 'Alignment Length', 
      'Mismatches', 'Gap Opens', 'Query Start', 'Query End', 
      'Subject Start', 'Subject End', 'e-Value', 'Bit Score'
    ]
    blast_results_df = pd.read_csv( # read in BLAST results
      f'{self.species_path}/blast_results.csv', names=result_columns
    )

    # extract the UniProt ID from the subject column
    blast_results_df['Subject'] = blast_results_df['Subject'].str.split('|').str[1]

    # take out "-#" portion of the subject UniProt ID because these are isoforms and 
    # won't be mapped properly to gene symbols
    blast_results_df['Subject'] = blast_results_df['Subject'].str.split('-').str[0]
    
    # map subject UniProt IDs to gene symbols
    blast_results_df['Subject Gene Symbol'] = blast_results_df['Subject'].map(self.uniprot_id_to_gene_symbol_map)
    
    # create a quality score based on % identity, alignment length, and query length
    blast_results_df['Query Length'] = blast_results_df['Query'].map(self.source_length_map)
    blast_results_df['Quality Score'] = blast_results_df['Percentage Identity'] * (blast_results_df['Alignment Length'] / blast_results_df['Query Length'])

    best_scores = blast_results_df.groupby('Query')['Quality Score'].transform(max) == blast_results_df['Quality Score']
    blast_results_df = blast_results_df[best_scores]

    for i, row in blast_results_df.iterrows():
      self.source_gene_assignment[row['Query']] = row['Subject Gene Symbol']
      self.source_protein_assignment[row['Query']] = row['Subject']

    # check if there are any query duplicates - update map within _pepmatch_tiebreak
    if blast_results_df.duplicated(subset=['Query']).any():
      self._pepmatch_tiebreak(blast_results_df)


  def _pepmatch_tiebreak(self, blast_results_df: pd.DataFrame) -> None:
    """First, get any source antigens that have ties for the best match. Then,
    using the epitopes associated with the source antigen, search them in 
    the selected proteome and find the gene that has the most matches.
    
    Args:
      blast_results_df: DataFrame of source antigen BLAST results.
    """
    from collections import Counter

    # get any source antigens that have ties for the best match
    source_antigens_with_ties = blast_results_df[blast_results_df.duplicated(subset=['Query'])]['Query'].unique()

    # preprocess with PEPMatch
    gp_proteome = f'{self.species_path}/gp_proteome.fasta' if os.path.exists(f'{self.species_path}/gp_proteome.fasta') else ''
    Preprocessor(
      proteome = f'{self.species_path}/proteome.fasta',
      preprocessed_files_path = f'{self.species_path}', 
      gene_priority_proteome = gp_proteome
    ).sql_proteome(k = 5)

    for source_antigen in source_antigens_with_ties:
      try:
        # get the epitopes associated with the source antigen
        epitopes = self.source_to_epitopes_map[source_antigen]
      except KeyError:
        # if there are no epitopes, then assign the gene and id to the first
        # blast match in blast_results_df of that source_antigen
        self.source_gene_assignment[source_antigen] = blast_results_df[blast_results_df['Query'] == source_antigen]['Subject Gene Symbol'].iloc[0]
        self.source_protein_assignment[source_antigen] = blast_results_df[blast_results_df['Query'] == source_antigen]['Subject'].iloc[0]
        continue

      matches_df = Matcher(
        query = epitopes, 
        proteome_file = f'{self.species_path}/proteome.fasta',
        max_mismatches = 0, 
        k = 5, 
        preprocessed_files_path = f'{self.species_path}',
        output_format='dataframe'
      ).match()

      if matches_df['Matched Sequence'].isna().all():
        # if there are no matches, then assign the gene and id to the first blast match
        self.source_gene_assignment[source_antigen] = blast_results_df[blast_results_df['Query'] == source_antigen]['Subject Gene Symbol'].iloc[0]
        self.source_protein_assignment[source_antigen] = blast_results_df[blast_results_df['Query'] == source_antigen]['Subject'].iloc[0]
        continue

      # get the uniprot id and gene symbol that has the most matches
      uniprot_id = Counter(matches_df['Protein ID']).most_common()[0][0]
      gene = Counter(matches_df['Gene']).most_common()[0][0]

      # update maps with the gene and id that has the most matches
      self.source_gene_assignment[source_antigen] = gene
      self.source_protein_assignment[source_antigen] = uniprot_id


  def _remove_files(self) -> None:
    """Delete all the files that were created when making the BLAST database."""
    for extension in ['pdb', 'phr', 'pin', 'pjs', 'pot', 'psq', 'ptf', 'pto']:
      try: # if DB wasn't create this will fail, so just pass
        os.remove(glob.glob(f'{self.species_path}/*.{extension}')[0])
      except IndexError:
        pass # the species will get all 0s for results
    
    # remove BLAST results and sources.fasta
    os.remove(f'{self.species_path}/blast_results.csv')
    os.remove(f'{self.species_path}/sources.fasta')

    # if proteome.db exists, remove it
    try:
      os.remove(f'{self.species_path}/proteome.db')
    except OSError:
      pass


  def _assign_allergens(self) -> dict:
    """Get allergen data from IUIS and create map."""
    url = 'http://www.allergen.org/csv.php?table=joint'
    allergen_df = pd.read_csv(url)
    allergen_map = allergen_df.set_index('AccProtein')['Name'].to_dict()

    for k, v in self.source_protein_assignment.items():
      if v in allergen_map.keys():
        self.source_protein_assignment[k] = allergen_map[v]


  def _get_manual_data(self) -> dict:
    """Get manual assignments from manual_assignments.csv and create map."""
    # manual_assignments.csv should be in the directory above this one
    directory = Path(__file__).resolve().parent.parent
    manual_df = pd.read_csv(directory / 'manual_assignments.csv')
    manual_map = dict(
      zip(
        manual_df['Accession'],
        manual_df['Parent Accession']
      )
    )
    return manual_map


def main():
  import argparse
  from get_data import DataFetcher

  parser = argparse.ArgumentParser()

  parser.add_argument(
    '-u', '--user',
    required=True,
    help='User for IEDB MySQL connection.')
  parser.add_argument(
    '-p', '--password',
    required=True,
    help='Password for IEDB MySQL connection.')
  parser.add_argument(
    '-a', '--all_species',
    action='store_true',
    help='Build protein tree for all IEDB species.')
  parser.add_argument(
    '-t', '--taxon_id',
    help='Taxon ID for the species to run protein tree.')
  
  args = parser.parse_args()
  user = args.user
  password = args.password
  all_species = args.all_species
  taxon_id = args.taxon_id

  species_df = pd.read_csv('species.csv')
  valid_taxon_ids = species_df['Species Taxon ID'].astype(str).tolist()

  # taxa, species name, and is_vertebrate mapppings
  all_taxa_map = dict(
    zip(
      species_df['Species Taxon ID'].astype(str),
      species_df['All Taxa']
    )
  )
  
  species_name_map = dict(
    zip(
      species_df['Species Taxon ID'].astype(str),
      species_df['Species Name']
    )
  )
  
  is_vertebrate_map = dict(
    zip(
      species_df['Species Taxon ID'].astype(str), 
      species_df['Is Vertebrate']
    )
  )

  if all_species:
    proteomes = {}
    for taxon_id in valid_taxon_ids:
      all_taxa = all_taxa_map[taxon_id]
      species_name = species_name_map[taxon_id]
      is_vertebrate = is_vertebrate_map[taxon_id]

      Fetcher = DataFetcher(user, password)
      epitopes_df = Fetcher.get_epitopes(all_taxa)
      sources_df = Fetcher.get_sources(all_taxa)

      if epitopes_df.empty or sources_df.empty:
        print(f'No epitopes for {species_name} ({taxon_id}). Skipping...\n')
        continue
      
      print(f'Assigning genes for {species_name} ({taxon_id})...')
      Assigner = GeneAndProteinAssigner(taxon_id, species_name, is_vertebrate)
      assigner_data = Assigner.assign(sources_df, epitopes_df)
      print('Done assigning genes.\n')

      print(f'Number of sources: {assigner_data[0]}')
      print(f'Number of epitopes: {assigner_data[1]}')
      print(f'Number of sources with BLAST matches: {assigner_data[2]}')
      print(f'Number of epitopes with a match: {assigner_data[3]}')
      print(f'Successful gene assignments: {(assigner_data[2] / assigner_data[0])*100:.1f}%')
      print(f'Successful parent assignments: {(assigner_data[3] / assigner_data[1])*100:.1f}%\n')

  else: # one species at a time
    assert taxon_id in valid_taxon_ids, f'{taxon_id} is not a valid taxon ID.'
    
    all_taxa = all_taxa_map[taxon_id]
    species_name = species_name_map[taxon_id]
    is_vertebrate = is_vertebrate_map[taxon_id]
    
    Fetcher = DataFetcher(user, password)
    epitopes_df = Fetcher.get_epitopes(all_taxa)
    sources_df = Fetcher.get_sources(all_taxa)

    assert not sources_df.empty, 'This species has no source antigens.'
    assert not epitopes_df.empty, 'This species has no epitopes.'

    print(f'Assigning genes for {species_name} ({taxon_id})...\n')
    Assigner = GeneAndProteinAssigner(taxon_id, species_name, is_vertebrate)
    assigner_data = Assigner.assign(sources_df, epitopes_df)
    print('Done assigning genes.\n')
    
    print(f'Number of sources: {assigner_data[0]}')
    print(f'Number of epitopes: {assigner_data[1]}')
    print(f'Number of sources with BLAST matches: {assigner_data[2]}')
    print(f'Number of epitopes with a match: {assigner_data[3]}')
    print(f'Successful gene assignments: {(assigner_data[2] / assigner_data[0])*100:.1f}%')
    print(f'Successful parent assignments: {(assigner_data[3] / assigner_data[1])*100:.1f}%\n')

if __name__ == '__main__':  
  main()