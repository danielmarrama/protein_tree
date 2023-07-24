#!/usr/bin/env python3

import warnings
warnings.filterwarnings('ignore')

import os
import glob
import multiprocessing
import pandas as pd

from ARC.classifier import SeqClassifier
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pathlib import Path
from pepmatch import Preprocessor, Matcher


class GeneAndProteinAssigner:
  def __init__(
    self,
    taxon_id,
    species_name,
    is_vertebrate,
    data_path = Path(__file__).parent.parent / 'data',
    bin_path = Path(__file__).parent.parent / 'bin',
    num_threads = multiprocessing.cpu_count() - 2
  ):
    self.species_dir = data_path / 'species' / f'{taxon_id}-{species_name.replace(" ", "_")}'
    self.taxon_id = taxon_id
    self.is_vertebrate = is_vertebrate
    self.data_path = data_path
    self.bin_path = bin_path
    self.num_threads = num_threads
    
    # initialize maps for assignments
    self.source_gene_assignment = {}
    self.source_protein_assignment = {}
    self.source_assignment_score = {}

    # create UniProt ID -> gene map
    self.proteome = pd.read_csv(f'{self.species_dir}/proteome.csv')
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
    num_matched_epitopes = self._assign_parents()
    self._assign_allergens()
    self._assign_manuals()

    # map source antigens to their best blast matches (UniProt ID and gene) for sources
    sources_df.loc[:, 'Assigned Gene'] = sources_df['Accession'].map(self.source_gene_assignment)
    sources_df.loc[:, 'Assigned Protein ID'] = sources_df['Accession'].map(self.source_protein_assignment)
    sources_df.loc[:, 'Assigned Protein Name'] = sources_df['Assigned Protein ID'].map(self.uniprot_id_to_name_map)
    sources_df.loc[:, 'Assignment Score'] = sources_df['Accession'].map(self.source_assignment_score)

    # map source antigens to their best blast matches (gene) for epitopes
    epitopes_df.loc[:, 'Assigned Gene'] = epitopes_df['Source Accession'].map(self.source_gene_assignment)
    epitopes_df.loc[:, 'Assigned Parent Protein ID'] = epitopes_df['Sequence'].map(self.epitope_parent_assignment)

    epitopes_df.drop_duplicates(subset='Sequence', inplace=True) # drop duplicate epitopes
    sources_df.drop(columns=['Sequence'], inplace=True) # drop sequence column for output
    
    self._remove_files()
    
    assigner_data = (
      num_sources,
      num_epitopes,
      num_matched_sources,
      num_matched_epitopes
    )

    return assigner_data, epitopes_df, sources_df


  def _assign_genes(
    self, sources_df: pd.DataFrame, epitopes_df: pd.DataFrame, num_sources: int
  ) -> None:
    """Assign a gene to the source antigens of a species.

    Run ARC for vertebrates to assign MHC/TCR/Ig to source antigens first.
    Run BLAST for all other source antigens to assign a gene and protein. 
    If there are ties, use PEPMatch to search the epitopes within the protein 
    sequences andselect the protein with the most matches.

    Args:
      sources_df: DataFrame of source antigens for a species.
      epitopes_df: DataFrame of epitopes for a species.
    """    
    self._sources_to_fasta(sources_df) # write sources to FASTA file

    if self.is_vertebrate:
      self._run_arc()

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
    with open(f'{self.species_dir}/sources.fasta', 'w') as f:
      SeqIO.write(seq_records, f, 'fasta')


  def _preprocess_proteome_if_needed(self) -> None:
    """Preprocess the proteome if the preprocessed files don't exist."""
    if not os.path.exists(f'{self.species_dir}/proteome.db'):
        gp_proteome = f'{self.species_dir}/gp_proteome.fasta' if os.path.exists(f'{self.species_dir}/gp_proteome.fasta') else ''
        Preprocessor(
          proteome = f'{self.species_dir}/proteome.fasta',
          preprocessed_files_path = f'{self.species_dir}',
          gene_priority_proteome=gp_proteome
        ).sql_proteome(k = 5)


  def _search_epitopes(self, epitopes: list, best_match: bool = True) -> pd.DataFrame:
    """Search epitopes within the proteome using PEPMatch."""
    df = Matcher(
      query = epitopes,
      proteome_file = f'{self.species_dir}/proteome.fasta',
      max_mismatches = 0, 
      k = 5, 
      preprocessed_files_path = f'{self.species_dir}', 
      best_match=best_match, 
      output_format='dataframe',
      sequence_version=False
    ).match()
    return df


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
    species_path = str(self.species_dir).replace('(', '\\(').replace(')', '\\)')

    os.system( # make BLAST database from proteome
      f'{self.bin_path}/makeblastdb -in {species_path}/proteome.fasta '\
      f'-dbtype prot > /dev/null'
    ) 
    os.system( # run blastp
      f'{self.bin_path}/blastp -query {species_path}/sources.fasta '\
      f'-db {species_path}/proteome.fasta '\
      f'-evalue 1 -num_threads {self.num_threads} -outfmt 10 '\
      f'-out {species_path}/blast_results.csv'
    ) 
    result_columns = [
      'Query', 'Target', 'Percentage Identity', 'Alignment Length', 
      'Mismatches', 'Gap Opens', 'Query Start', 'Query End', 
      'Target Start', 'Target End', 'e-Value', 'Bit Score'
    ]
    blast_results_df = pd.read_csv( # read in BLAST results
      f'{self.species_dir}/blast_results.csv', names=result_columns
    )

    # extract the UniProt ID from the target column
    blast_results_df['Target'] = blast_results_df['Target'].str.split('|').str[1]

    # take out "-#" portion of the target UniProt ID because these are isoforms and 
    # won't be mapped properly to gene symbols
    blast_results_df['Target'] = blast_results_df['Target'].str.split('-').str[0]
    
    # map target UniProt IDs to gene symbols
    blast_results_df['Target Gene Symbol'] = blast_results_df['Target'].map(self.uniprot_id_to_gene_map)
    
    # create a quality score based on % identity, alignment length, and query length
    blast_results_df['Query Length'] = blast_results_df['Query'].map(self.source_length_map)
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
      self.source_gene_assignment[row['Query']] = row['Target Gene Symbol']
      self.source_protein_assignment[row['Query']] = row['Target']
      self.source_assignment_score[row['Query']] = row['Quality Score']


  def _assign_parents(self) -> tuple:
    """Assign a parent protein to each epitope.
    
    Preprocess the proteome and then search all the epitopes within
    the proteome using PEPMatch. Then, assign the parent protein
    to each epitope by selecting the best isoform of the assigned gene for
    its source antigen.
    """
    self._preprocess_proteome_if_needed()

    all_epitopes = []
    for epitopes in self.source_to_epitopes_map.values():
      all_epitopes.extend(epitopes)

    all_matches_df = self._search_epitopes(all_epitopes, best_match=False)

    self.epitope_parent_assignment = {}
    for antigen, epitopes in self.source_to_epitopes_map.items():
      for epitope in epitopes:

        # check if the epitope is found in the protein the antigen is assigned to
        matches_df = all_matches_df[all_matches_df['Query Sequence'] == epitope]
        try:
          epitope_protein = matches_df[matches_df['Protein ID'] == self.source_protein_assignment[antigen]]
        except KeyError: # if the source antigen has no assignment, skip
          continue
      
        if epitope_protein.empty: # otherwise, get an isoform of the assigned gene
          epitope_protein = matches_df[matches_df['Gene'] == self.source_gene_assignment[antigen]]

          if len(epitope_protein) > 1: # if there are ties, select the protein with the best protein existence level
            epitope_protein = epitope_protein.loc[epitope_protein['Protein Existence Level'].idxmin()]

        if epitope_protein.empty: # still empty, skip
          continue

        if isinstance(epitope_protein['Protein ID'], pd.Series):
            self.epitope_parent_assignment[epitope] = epitope_protein['Protein ID'].values[0]
        else: # sometimes epitope_protein['Protein ID'] is a string
            self.epitope_parent_assignment[epitope] = epitope_protein['Protein ID']
    
    # count the number of epitopes with matches
    num_matched_epitopes = len(self.epitope_parent_assignment.keys())
    
    return num_matched_epitopes


  def _run_arc(self) -> None:
    """Run ARC to assign MHC/TCR/Ig to source antigens."""

    print('Running ARC for MHC/TCR/Ig assignments...\n')

    # pass sources.fasta to ARC
    SeqClassifier(
      outfile = f'{self.species_dir}/ARC_results.tsv',
      threads = self.num_threads,
      hmmer_path = str(self.bin_path) + '/',
      blast_path = str(self.bin_path) + '/',
    ).classify_seqfile(f'{self.species_dir}/sources.fasta')

    arc_results = pd.read_csv(f'{self.species_dir}/ARC_results.tsv', sep='\t')
    
    if not arc_results.dropna(subset=['class']).empty:
      arc_assignments = arc_results.set_index('id')['class'].to_dict()
      self.source_gene_assignment.update(arc_assignments)


  def _assign_allergens(self) -> None:
    """Get allergen data from allergen.org and then assign allergens to sources."""
    allergen_df = pd.read_csv(self.data_path / 'allergens.csv')
    allergen_map = allergen_df.set_index('AccProtein')['Name'].to_dict()

    for k, v in self.source_protein_assignment.items():
      if v in allergen_map.keys():
        self.uniprot_id_to_name_map[v] = allergen_map[v]


  def _assign_manuals(self) -> None:
    """Get manual assignments from manual_assignments.csv and then assign
    genes and proteins to sources.
    """
    # manual_assignments.csv should be in the directory above this one
    manual_df = pd.read_csv(self.data_path / 'manual_assignments.csv')
    manual_gene_map = manual_df.set_index('Accession')['Accession Gene'].to_dict()
    manual_protein_id_map = manual_df.set_index('Accession')['Parent Accession'].to_dict()
    manual_protein_name_map = manual_df.set_index('Accession')['Parent Name'].to_dict()

    for k, v in self.source_gene_assignment.items():
      if k in manual_gene_map.keys():
        self.source_gene_assignment[k] = manual_gene_map[k]
        self.source_protein_assignment[k] = manual_protein_id_map[k]
        self.uniprot_id_to_name_map[k] = manual_protein_name_map[k]


  def _remove_files(self) -> None:
    """Delete all the files that were created when making the BLAST database."""
    for extension in ['pdb', 'phr', 'pin', 'pjs', 'pot', 'psq', 'ptf', 'pto']:
      try: # if DB wasn't created this will throw an error
        os.remove(glob.glob(f'{self.species_dir}/*.{extension}')[0])
      except IndexError:
        pass 
    
    os.remove(f'{self.species_dir}/blast_results.csv')
    os.remove(f'{self.species_dir}/sources.fasta')

    if self.is_vertebrate:
      os.remove(f'{self.species_dir}/ARC_results.tsv')

    # if proteome.db exists, remove it
    try:
      os.remove(f'{self.species_dir}/proteome.db')
    except OSError:
      pass

