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

    # create species path with species taxon and name; example: "24-Shewanella putrefaciens"
    species_id_to_name_map = dict(zip(self.species_df['Taxon ID'].astype(str), self.species_df['Species Label']))
    self.species_path = f'species/{taxon_id}-{species_id_to_name_map[taxon_id].replace(" ", "_")}'

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

    # create blast database
    self._create_blast_db()

    # run blast
    self._run_blast()

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
      sources_df[sources_df['Sequence'].isna()].to_csv(f'{self.species_path}/sources_missing_seqs.csv', index=False)
   
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

    with open(f'{self.species_path}/sources.fasta', 'w') as f:
      SeqIO.write(seq_records, f, 'fasta')

  def _create_blast_db(self):
    os.system(f'./makeblastdb -in {self.species_path}/proteome.fasta -dbtype prot')

  def _run_blast(self):
    '''
    BLAST source antigens against the selected proteome, then read in with
    pandas and assign column names. By default, blastp doesn't return header.
    '''
    os.system(f'./blastp -query {self.species_path}/sources.fasta '\
              f'-db {self.species_path}/proteome.fasta '\
              f'-evalue 0.0001 -num_threads 1 -outfmt 10 '\
              f'-out {self.species_path}/blast_results.csv'
    )
    
    result_columns = ['query', 'subject', '%_identity', 'alignment_length', 
                      'mismatches', 'gap_opens', 'q_start', 'q_end', 
                      's_start', 's_end', 'evalue', 'bit_score']

    # read in results that were just written and rewrite with column headers
    pd.read_csv(f'{self.species_path}/blast_results.csv', names=result_columns).to_csv(f'{self.species_path}/blast_results.csv', index=False)


  def _no_blast_match(self):
    '''Write sources that have no BLAST match to a file.'''

    # get all source ids
    source_ids = []
    for record in list(SeqIO.parse(f'{self.species_path}/sources.fasta', 'fasta')):
      source_ids.append(str(record.id))
    
    # get BLAST results and then get ids that are not in results
    blast_results = pd.read_csv(f'{self.species_path}/blast_results.csv')
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