#!/usr/bin/env python3

import argparse
from select_proteome import ProteomeSelector
from assign_genes import GeneAssigner

# definite command line args which will take in a taxon ID
# and run the select_proteome.py script as well as  the 
# assign_gene.py script to get the gene assignments and parent proteins
# in order to build the protein tree for each species
def main():
  parser = argparse.ArgumentParser(description='Build protein tree for a species.')
  parser.add_argument('taxon_id', type=int, help='Taxon ID of species.')
  args = parser.parse_args()

  # run select_proteome.py script


  # run assign_gene.py script

if __name__ == '__main__':
  main()