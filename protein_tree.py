#!/usr/bin/env python3

import os
import argparse
import pandas as pd 

import get_data
import select_proteome
import assign_genes


def run_protein_tree(taxon_id):
  """
  Build protein tree for a species.
  """
  # make species folder if it doesn't exist
  os.makedirs(f'species/{taxon_id}', exist_ok=True)

  print('Getting data from IEDB MySQL backend...')
  epitopes = get_data.get_epitopes(taxon_id)
  sources = get_data.get_sources(taxon_id)
  print('Done getting data.\n')

  print('Selecting the best proteome...')
  proteome_id, proteome_taxon, proteome_type = select_proteome.select_proteome(taxon_id)
  print(f'Got the best proteome.\n')

  # print('Assigning genes to proteins...')
  # assign_genes.assign_genes(taxon_id, proteome_id)

def main():
  parser = argparse.ArgumentParser(description='Build protein tree for a species.')
  parser.add_argument('taxon_id', help='Taxon ID of species.')
  args = parser.parse_args()

  taxon_id = args.taxon_id

  # TODO: replace species.csv with a call to the MySQL backend
  species_df = pd.read_csv('species.csv')
  valid_taxon_ids = species_df['Taxon ID'].astype(str).tolist()
  id_to_names = dict(zip(species_df['Taxon ID'].astype(str), species_df['Species Label']))

  # make species folder if it doesn't exist
  os.makedirs('species', exist_ok=True)

  if taxon_id == 'all':
    for taxon_id in valid_taxon_ids:
      print(f'Building protein tree for {id_to_names[taxon_id]} (ID: {taxon_id})...\n')
      run_protein_tree(taxon_id)
      print('Protein tree build done.')

  else:
    assert taxon_id in valid_taxon_ids, f'{taxon_id} is not a valid taxon ID.'
    print(f'Building protein tree for {id_to_names[taxon_id]} (ID: {taxon_id})...\n')
    run_protein_tree(taxon_id)
    print('Protein tree build done.')

if __name__ == '__main__':
  main()