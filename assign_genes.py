#!/usr/bin/env python3

import os
import glob
import pandas as pd

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

class GeneAssigner:
  def __init__(self, taxon_id):
    self.species_df = pd.read_csv('species.csv')
    self.taxon_id = taxon_id

    # create species path with species taxon and name; example: "24-Shewanella_putrefaciens"
    species_id_to_name_map = dict(zip(self.species_df['Taxon ID'].astype(str), self.species_df['Species Label']))
    self.species_path = f'species/{taxon_id}-{species_id_to_name_map[taxon_id].replace(" ", "_")}'

    # create UniProt ID to gene symbol map from proteome.csv file
    proteome = pd.read_csv(f'{self.species_path}/proteome.csv')
    self.uniprot_id_to_gene_symbol_map = dict(zip(proteome['UniProt ID'], proteome['Gene Symbol']))

  def assign_genes(self, sources_df, epitopes_df):
    """
    Assign a gene to the source antigens of a species.

    First, make a BLAST database of the selected proteome.
    Then, run blastp with the source antigens against the proteome.
    Based on the BLAST results, select the best protein match by
    sequence identity and the highest alignment length.

    If there are ties, use PEPMatch to search the epitopes within
    the protein sequences and select the protein with the most matches.
    """
    if sources_df.empty:
      return

    # create source to epitope map and write sources to FASTA file
    self.source_to_epitopes_map = self._create_source_to_epitopes_map(epitopes_df)
    num_sources, num_sources_missing_seqs = self._sources_to_fasta(sources_df)

    # create BLAST database, run blastp, and remove db files
    self._create_blast_db()
    blast_results_df = self._run_blast()

    # remove sources that don't have any BLAST matches and get the counts
    num_no_blast_matches, num_with_blast_matches = self._no_blast_matches()

    # get best blast matches for each source antigen
    self._get_best_blast_matches(blast_results_df)

    # map source antigens to their best blast matches (UniProt ID and gene)
    sources_df['Assigned Gene'] = sources_df['Accession'].map(self.best_blatch_match_gene_map)
    sources_df['Assigned UniProt ID'] = sources_df['Accession'].map(self.best_blatch_match_id_map)
    
    # write sources with assigned genes to file
    sources_df.to_csv(f'{self.species_path}/sources.csv', index=False)

    # remove blast DB and result files
    self._remove_files()

    return num_sources, num_sources_missing_seqs, num_no_blast_matches, num_with_blast_matches

  def assign_parents(self, epitopes_df):
    pass

  def _create_source_to_epitopes_map(self, epitopes_df):
    source_to_epitopes_map = {}
    for i, row in epitopes_df.iterrows():
      if row['Source Accession'] in source_to_epitopes_map.keys():
        source_to_epitopes_map[row['Source Accession']].append(row['Peptide'])
      else:
        source_to_epitopes_map[row['Source Accession']] = [row['Peptide']]
    
    return source_to_epitopes_map 

  def _sources_to_fasta(self, sources_df):
    """
    Write source antigens to FASTA file. If a source antigen is missing
    a sequence, write it to a separate file for logging.
    """  
    # write sources that are missing sequences to file and then drop those
    if not sources_df[sources_df['Sequence'].isna()].empty:
      num_sources_missing_seqs = len(sources_df[sources_df['Sequence'].isna()])
      sources_df[sources_df['Sequence'].isna()].to_csv(f'{self.species_path}/sources_missing_seqs.csv', index=False)
    else:
      num_sources_missing_seqs = 0

    sources_df.dropna(subset=['Sequence'], inplace=True)
        
    # create seq records of sources with ID and sequence
    # TEST: use 1,000 sources for testing
    # TODO: REMOVE THIS 
    seq_records = []
    for i, row in sources_df.iloc[0:1000, :].iterrows():
      seq_records.append(
        SeqRecord(
          Seq(row['Sequence']),
          id=row['Accession'],
          description='')
      )

    with open(f'{self.species_path}/sources.fasta', 'w') as f:
      SeqIO.write(seq_records, f, 'fasta')

    return len(sources_df), num_sources_missing_seqs

  def _create_blast_db(self):
    # escape parentheses in species path
    species_path = self.species_path.replace('(', '\(').replace(')', '\)')
    os.system(f'./makeblastdb -in {species_path}/proteome.fasta -dbtype prot')

  def _run_blast(self):
    '''
    BLAST source antigens against the selected proteome, then read in with
    pandas and assign column names. By default, blastp doesn't return header.
    '''
    # escape parentheses in species path
    species_path = self.species_path.replace('(', '\(').replace(')', '\)')  
    os.system(f'./blastp -query {species_path}/sources.fasta '\
              f'-db {species_path}/proteome.fasta '\
              f'-evalue 0.0001 -num_threads 1 -outfmt 10 '\
              f'-out {species_path}/blast_results.csv'
    )
    
    result_columns = ['Query', 'Subject', 'Percentage Identity', 'Alignment Length', 
                      'Mismatches', 'Gap Opens', 'Query Start', 'Query End', 
                      'Subject Start', 'Subject End', 'e-Value', 'Bit Score']

    # read in results that were just written and map subject UniProt IDs to gene symbols
    blast_results_df = pd.read_csv(f'{self.species_path}/blast_results.csv', names=result_columns)
    blast_results_df['Subject'] = blast_results_df['Subject'].str.split('|').str[1]
    blast_results_df['Subject Gene Symbol'] = blast_results_df['Subject'].map(self.uniprot_id_to_gene_symbol_map)
    
    # write results with column header and gene symbols to file
    blast_results_df.to_csv(f'{self.species_path}/blast_results.csv', index=False)

    return blast_results_df

  def _remove_files(self):
    """Delete all the files that were created when making the BLAST database."""
    for extension in ['pdb', 'phr', 'pin', 'pjs', 'pot', 'psq', 'ptf', 'pto']:
      os.remove(glob.glob(f'{self.species_path}/*.{extension}')[0])
    
    # remove BLAST results and sources.fasta
    os.remove(f'{self.species_path}/blast_results.csv')
    os.remove(f'{self.species_path}/sources.fasta')

    # if proteome.db exists, remove it
    try:
      os.remove(f'{self.species_path}/proteome.db')
    except OSError:
      pass

  def _no_blast_matches(self):
    '''Write sources that have no BLAST match to a file.'''
    # get all source antigen ids
    source_ids = []
    for record in list(SeqIO.parse(f'{self.species_path}/sources.fasta', 'fasta')):
      source_ids.append(str(record.id))
    
    # get BLAST results and then get ids that are not in results
    blast_results = pd.read_csv(f'{self.species_path}/blast_results.csv')
    blast_result_ids = list(blast_results['Query'].unique())

    no_blast_match_ids = list(set(source_ids) - set(blast_result_ids))

    # write no BLAST match ids to a file if there are any 
    if len(no_blast_match_ids) > 0:
      with open(f'{self.species_path}/no_blast_match_ids.txt', 'w') as f:
        for id in no_blast_match_ids:
          f.write(f'{id}\n')

    # return the number of sources that have BLAST matches and no BLAST matches
    return len(no_blast_match_ids), len(blast_result_ids)

  def _get_best_blast_matches(self, blast_results_df):
    """
    Get the best BLAST match for each source antigen by sequence identity and 
    alignment length. If there are multiple matches with the same % identity 
    and alignment length, then use _pepmatch_tiebreak to determine the best match.
    """
    # get the best match for each source antigen by % identity
    index = blast_results_df.groupby(['Query'])['Percentage Identity'].transform(max) == blast_results_df['Percentage Identity']
    blast_results_df = blast_results_df[index]

    # and alignment length
    index = blast_results_df.groupby(['Query'])['Alignment Length'].transform(max) == blast_results_df['Alignment Length']
    blast_results_df = blast_results_df[index]

    self.best_blatch_match_gene_map, self.best_blatch_match_id_map = {}, {}
    for i, row in blast_results_df.iterrows():
      self.best_blatch_match_gene_map[row['Query']] = row['Subject Gene Symbol']
      self.best_blatch_match_id_map[row['Query']] = row['Subject']

    # check if there are any query duplicates - update map within _pepmatch_tiebreak
    if blast_results_df.duplicated(subset=['Query']).any():
      self._pepmatch_tiebreak(blast_results_df)

  def _pepmatch_tiebreak(self, blast_results_df):
    """
    First, get any source antigens that have ties for the best match. Then,
    using the epitopes associated with the source antigen, search them in 
    the selected proteome and find the gene that has the most matches.
    """
    from pepmatch import Preprocessor, Matcher
    from collections import Counter

    # get any source antigens that have ties for the best match
    source_antigens_with_ties = blast_results_df[blast_results_df.duplicated(subset=['Query'])]['Query'].unique()

    for source_antigen in source_antigens_with_ties:
      # get the epitopes associated with the source antigen
      epitopes = self.source_to_epitopes_map[source_antigen]

      # search the epitopes in the selected proteome
      Preprocessor(f'{self.species_path}/proteome.fasta', 'sql', f'{self.species_path}').preprocess(k=5)
      matches_df = Matcher(epitopes, 'proteome', 0, 5, f'{self.species_path}', output_format='dataframe').match()
      
      # get the uniprot id and gene symbol that has the most matches
      uniprot_id = Counter(matches_df['Protein ID']).most_common()[0][0]
      gene = Counter(matches_df['Gene']).most_common()[0][0]

      # update maps with the gene and id that has the most matches
      self.best_blatch_match_gene_map[source_antigen] = gene
      self.best_blatch_match_id_map[source_antigen] = uniprot_id

def main():
  import argparse
  from get_data import DataFetcher

  parser = argparse.ArgumentParser()

  parser.add_argument('-u', '--user', required=True, help='User for IEDB MySQL connection.')
  parser.add_argument('-p', '--password', required=True, help='Password for IEDB MySQL connection.')
  parser.add_argument('-a', '--all_species', action='store_true', help='Build protein tree for all IEDB species.')
  parser.add_argument('-t', '--taxon_id', help='Taxon ID for the species to run protein tree.')
  
  args = parser.parse_args()
  user = args.user
  password = args.password
  all_species = args.all_species
  taxon_id = args.taxon_id

  # read in IEDB species data
  species_df = pd.read_csv('species.csv')
  valid_taxon_ids = species_df['Taxon ID'].astype(str).tolist()

  # dicts for mapping taxon IDs to all their taxa and their names
  all_taxa_map = dict(zip(species_df['Taxon ID'].astype(str), species_df['All Taxa']))
  species_id_to_name_map = dict(zip(species_df['Taxon ID'].astype(str), species_df['Species Label']))

  # do proteome selection for all IEDB species
  if taxon_id == 'all':
    proteomes = {}
    for taxon_id in valid_taxon_ids:
      # get data for taxon ID
      Fetcher = DataFetcher(user, password, taxon_id, all_taxa_map[taxon_id])
      epitopes_df = Fetcher.get_epitopes()
      sources_df = Fetcher.get_sources()
      
      Assigner = GeneAssigner(taxon_id)
      num_sources, num_sources_missing_seqs, num_no_blast_matches, num_with_blast_matches = Assigner.assign_genes(sources_df, epitopes_df)
      # Assigner.assign_parents()

  # or just one species at a time - check if its valid
  else:
    assert taxon_id in valid_taxon_ids, f'{taxon_id} is not a valid taxon ID.'
    # get data for taxon ID
    Fetcher = DataFetcher(user, password, taxon_id, all_taxa_map[taxon_id])
    epitopes_df = Fetcher.get_epitopes()
    sources_df = Fetcher.get_sources()

    assert not sources_df.empty, 'This species has no source antigens.'

    Assigner = GeneAssigner(taxon_id)
    num_sources, num_sources_missing_seqs, num_no_blast_matches, num_with_blast_matches = Assigner.assign_genes(sources_df, epitopes_df)
    # Assigner.assign_parents()

if __name__ == '__main__':  
  main()