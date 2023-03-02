#!/usr/bin/env python3

import os

import pandas as pd

from pepmatch import Preprocessor, Matcher

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def sources_to_fasta(taxon_id):
  # sources that are missing sequences - write to file and then drop those
  sources[self.sources['Source Sequence'].isna()].to_csv(f'species/{taxon_id}/sources_missing_seqs.csv', index=False)
  sources.dropna(subset=['Source Sequence'], inplace=True)

  # create seq records of sources with ID and sequence
  # TEST: use 1,000 sources for testing
  seq_records = []
  for i, row in self.sources.iloc[0:1000, :].iterrows():
    seq_records.append(
      SeqRecord(
        row['Source Sequence'], 
        id=row['Source Accession'], 
        description='')
    )

  with open(f'species/{self.taxon_id}/sources.fasta', 'w') as f:
    SeqIO.write(seq_records, f, 'fasta')

def create_source_epitope_map(self):
  source_epitope_map = {}
  for i, row in self.epitopes.iterrows():
    if row['Source Accession'] in source_epitope_map.keys():
      source_epitope_map[row['Source Accession']].append(row['Peptide'])
    else:
      source_epitope_map[row['Source Accession']] = [row['Peptide']]
  
  return source_epitope_map 

def create_blast_db(self):
  os.system(f'./makeblastdb -in {self.taxon_id}/proteome.fasta -dbtype prot')

def run_blast(self):
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


def no_blast_match(self):
  '''Write sources that have no BLAST match to a file.'''

  # get all source ids
  source_ids = []
  for record in list(SeqIO.parse(f'{self.taxon_id}/sources.fasta', 'fasta')):
    source_ids.append(str(record.id))
  
  # get BLAST results and then get ids that are not in results
  blast_results = pd.read_csv(f'{self.taxon_id}/blast_results.csv')
  blast_result_ids = list(blast_results['query'].unique())

def pepmatch_tiebreak():
  pass

def assign_genes(taxon_id):
  '''
  Assign gene to sources based on BLAST results.
  '''

  # read in proteome and create a map of id to gene
  proteome = pd.read_csv(f'{self.taxon_id}/proteome.csv')
  id_to_gene_map = dict(zip(proteome['id'], hproteomep['gene']))

  # read in BLAST results
  df = pd.read_csv(f'{self.taxon_id}/blast_results.csv')

  # df['subject_gene'] = df['subject'].str.split('|').str[1].map(gene_id_map)

  # df.drop_duplicates(subset=['query', 'subject_gene'], inplace=True)


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