#!/usr/bin/env python3

import os

import pandas as pd

from pepmatch import Preprocessor, Matcher

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


class GeneAssignment(object):
  pass

def pepmatch_tiebreak():
  pass

def assign_genes():
  pass

def run_blast(taxon_id):
  '''
  Run blastp with source antigens to the proteome. 
  '''
  os.system('./blastp -query test.fasta -db %s/proteome.fasta'\
            '-evalue 0.0001 -num_threads 14 -outfmt 10 -out'\
            '%s/blast_results.csv' % (taxon_id, taxon_id)
  )

def create_blast_db(taxon_id):
  os.system('makeblastdb -in %s/proteome.fasta -dbtype prot' % taxon_id)

def create_source_epitope_map(sources, epitopes):
  source_epitope_map = {}
  for i, row in epitopes.iterrows():
    if row['Source Accession'] in source_epitope_map.keys():
      source_epitope_map[row['Source Accession']].append(row['Peptide'])
    else:
      source_epitope_map[row['Source Accession']] = [row['Peptide']]
  
  return source_epitope_map 

def sources_to_fasta(sources):
  # sources that are missing sequences - write to file and then drop those
  sources[sources['Source Sequence'].isna()].to_csv('%s/sources_missing_seqs.csv', index=False)
  sources.dropna(subset=['Source Sequence'], inplace=True)

  # create seq records of sources with ID and sequence
  # TEST: use 10,000 sources for testing
  seq_records = []
  for i, row in sources.iloc[0:10000, :].iterrows():
    seq_records.append(
      SeqRecord(
        row['Source Sequence'], 
        id=row['Source Accession'], 
        description='')
    )

  with open('%s/sources.fasta', 'w') as f:
    SeqIO.write(seq_records, f, 'fasta')

def get_data(taxon_id):
    sources = pd.read_csv('%s/sources.tsv' % taxon_id, sep='\t')
    epitopes =  pd.read_csv('%s/epitopes.tsv' % taxon_id, sep='\t')
    return sources, epitopes


def run(taxon_id):
  # TODO: 
  # bucket sources that did not get a BLAST match
  # limit results to highest % identity / alignment length
  # tiebreak the tied results with PEPMatch epitope search

  # get sources and epitopes data
  sources, epitopes = get_data(taxon_id) 

  # write source antigens to FASTA with their IDs
  sources_to_fasta(sources)

  # create dict of source antigens to their list of epitopes
  source_epitope_map = create_source_epitope_map(sources, epitopes)

  # create a BLAST database of the species proteome
  create_blast_db(taxon_id)

if __name__ == '__main__':
  taxon_id = '9606'

  print('Running gene assignments on %s...\n' % taxon_id)

  run(taxon_id)