#!/usr/bin/env python3

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

def run_blast():
    pass

def create_blast_db():
    pass

def create_source_epitope_map(sources, epitopes):
    pass

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
    # write sources to FASTA file
    # create source -> epitope mapping
    # create BLAST DB on proteome 
    # run BLAST using sources on proteome
    # bucket sources that did not get a BLAST match
    # limit results to highest % identity / alignment length
    # tiebreak the tied results with PEPMatch epitope search

    sources, epitopes = get_data(taxon_id) 
    sources_to_fasta(sources)
    # source_epitope_map = create_source_epitope_map(sources, epitopes)

if __name__ == '__main__':
    taxon_id = '9606'
    run(taxon_id)