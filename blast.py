#!/usr/bin/env python3

import os, glob

import pandas as pd
import matplotlib.pyplot as plt

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Blast.Applications import NcbiblastpCommandline

from sklearn import metrics
from sklearn.metrics import roc_curve
from sklearn.metrics import roc_auc_score



def run_blast():
  blastx_cline = NcbiblastpCommandline(
              num_threads=14,
              cmd='./blastp',
              query = 'sources.fasta', # 10,000 antigens as test
              db = '9606_reference_with_isoforms.fasta', 
              evalue=1.0E-5, outfmt=10, out='blast_output.csv')

  stdout, stderr = blastx_cline()

  # remove_files(fold)

if __name__ == '__main__':
  run_blast()






















  # def remove_files(fold):
#   '''
#   Removes db.fasta file and BLAST DB files created for each fold.
#   '''
#   os.remove('db.fasta')
#   os.remove('fold%i.fasta' % fold)
#   for extension in ['pdb', 'phr', 'pin', 'psq', 'ptf', 'pot', 'pto', 'pjs']:
#     os.remove(glob.glob('*.' + extension)[0])
