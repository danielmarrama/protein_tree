#!/usr/bin/env python3

import os
import argparse
import pandas as pd 

import select_proteome

def main():
  parser = argparse.ArgumentParser(description='Build protein tree for a species.')
  parser.add_argument('taxon_id', help='Taxon ID of species.')
  args = parser.parse_args()

  taxon_id = args.taxon_id

  species_df = pd.read_csv('species.csv')
  valid_taxon_ids = species_df['Taxon ID'].astype(str).tolist()

  if not taxon_id == 'all':
    assert taxon_id in valid_taxon_ids, f'{taxon_id} is not a valid taxon ID.'
  else:
    raise ValueError('Taxon ID is not valid.')

  # make species folder if it doesn't exist
  os.makedirs('species', exist_ok=True)

  # run select_proteome.py script
  print(select_proteome.select_proteome(taxon_id))

  # run assign_gene.py script

if __name__ == '__main__':
  main()