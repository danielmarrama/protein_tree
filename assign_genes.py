#!/usr/bin/env python3

import os

import pandas as pd

from pepmatch import Preprocessor, Matcher

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


class GeneAssigner:
  def __init__(self, taxon_id):
    self.species_df = pd.read_csv('species.csv')
    self.taxon_id = taxon_id

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
    # create source to epitope map
    self.source_to_epitopes_map = self._create_source_to_epitopes_map(epitopes_df)

    # write sources to fasta
    self._sources_to_fasta(sources_df)

    # # create blast database
    # self._create_blast_db()

    # # run blast
    # self.run_blast()

    # # read in blast results
    # self.blast_results = self.read_blast_results()

    # # assign genes to sources
    # self.assign_genes_to_sources()

    # # assign parents to sources
    # self.assign_parents()

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
      sources_df[sources_df['Sequence'].isna()].to_csv(f'species/{self.taxon_id}/sources_missing_seqs.csv', index=False)
   
    sources_df.dropna(subset=['Sequence'], inplace=True)

    # create seq records of sources with ID and sequence
    # TEST: use 1,000 sources for testing
    seq_records = []
    for i, row in sources_df.iloc[0:1000, :].iterrows():
      seq_records.append(
        SeqRecord(
          Seq(row['Sequence']),
          id=row['Accession'],
          description='')
      )

    with open(f'species/{self.taxon_id}/sources.fasta', 'w') as f:
      SeqIO.write(seq_records, f, 'fasta')

  def _create_blast_db(self):
    os.system(f'./makeblastdb -in species/{self.taxon_id}/proteome.fasta -dbtype prot')

  def _run_blast(self):
    '''
    Run blastp with source antigens to the proteome.
    Then read in with pandas and assign column names.
    '''
    path = f'{self.taxon_id}/blast_results.csv'
    os.system(f'./blastp -query test.fasta -db {self.taxon_id}/proteome.fasta'\
              f' -evalue 0.0001 -num_threads 14 -outfmt 10 -out {path}'
    )
    
    blast_columns = ['query', 'subject', '%_identity', 'alignment_length', 
                      'mismatches', 'gap_opens', 'q_start', 'q_end', 
                      's_start', 's_end', 'evalue', 'bit_score']
    pd.read_csv(path, names=blast_columns).to_csv(path, index=False)


  def _no_blast_match(self):
    '''Write sources that have no BLAST match to a file.'''

    # get all source ids
    source_ids = []
    for record in list(SeqIO.parse(f'{self.taxon_id}/sources.fasta', 'fasta')):
      source_ids.append(str(record.id))
    
    # get BLAST results and then get ids that are not in results
    blast_results = pd.read_csv(f'{self.taxon_id}/blast_results.csv')
    blast_result_ids = list(blast_results['query'].unique())

  def _pepmatch_tiebreak():
    pass


def main():
  # define command line args which will take in a taxon ID
  parser = argparse.ArgumentParser(description='Select the taxon ID to assign the source antigen to genes.')
  parser.add_argument('taxon_id', help='Taxon ID for the species.')
  args = parser.parse_args()
  taxon_id = args.taxon_id

  # TODO: 
  # limit results to highest % identity / alignment length
  # tiebreak the tied results with PEPMatch epitope search

  # write source antigens to FASTA with their IDs
  print('Writing sources to FASTA file...')
  sources_to_fasta(taxon_id)
  print('Done.')

  # create dict of source antigens to their list of epitopes
  print('Creating source to epitopes mapping...')
  source_epitope_map = self.create_source_epitope_map()
  print('Done.')

  # create a BLAST database of the species proteome
  print('Creating BLAST database of proteome...')
  self.create_blast_db()
  print('Done.')
  
  # BLAST sources against proteome  
  print('BLASTing sources against proteome database...')
  self.run_blast()
  print('Done.')

  # bucket sources that did not get a BLAST match
  print('Writing sources that did not BLAST match to file...')
  self.no_blast_match()
  print('Done.')

if __name__ == '__main__':  
  main()