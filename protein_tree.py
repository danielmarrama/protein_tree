#!/usr/bin/env python3

import os
import argparse
import pandas as pd 

from get_data import DataFetcher
import select_proteome
import assign_genes


def run_protein_tree(user, password, taxon_id):
  """
  Build protein tree for a species.
  """
  # make species folder if it doesn't exist
  os.makedirs(f'species/{taxon_id}', exist_ok=True)

  print('Getting epitopes and sources data...')
  fetcher = DataFetcher(user, password, taxon_id)
  epitopes_df = fetcher.get_epitopes()
  sources_df = fetcher.get_sources()
  print(epitopes_df)
  print('Done getting data.\n')

  print('Selecting the best proteome...\n')
  num_of_proteomes, proteome_id, proteome_taxon, proteome_type = select_proteome.select_proteome(taxon_id, epitopes_df)
  
  print(f'Number of proteomes: {num_of_proteomes}')
  print('Got the best proteome:')
  print(f'Proteome ID: {proteome_id}')
  print(f'Proteome taxon: {proteome_taxon}')
  print(f'Proteome type: {proteome_type}')

  print('\nAssigning genes to proteins...\n')
  # assign_genes.assign_genes(taxon_id, proteome_id)

def main():
  # define command line args which will take in a taxon ID
  parser = argparse.ArgumentParser()
  
  parser.add_argument('-u', '--user', required=True, help='User for IEDB MySQL connection.')
  parser.add_argument('-p', '--password', required=True, help='User for IEDB MySQL connection.')
  parser.add_argument('-t', '--taxon_id', required=True, help='Taxon ID for the species to pull data for.')
  
  args = parser.parse_args()
  user = args.user
  password = args.password
  taxon_id = args.taxon_id

  # TODO: replace species.csv with a call to the MySQL backend
  species_df = DataFetcher(user, password, taxon_id).get_species()
  valid_taxon_ids = species_df['Taxon ID'].astype(str).tolist()
  id_to_names = dict(zip(species_df['Taxon ID'].astype(str), species_df['Species Label']))

  # make species folder if it doesn't exist
  os.makedirs('species', exist_ok=True)

  # run protein tree for all IEDB species 
  if taxon_id == 'all':
    for taxon_id in valid_taxon_ids:
      print(f'Building protein tree for {id_to_names[taxon_id]} (ID: {taxon_id})...\n')
      run_protein_tree(user, password, taxon_id)
      print('Protein tree build done.')

  # or one species at a time
  else:
    assert taxon_id in valid_taxon_ids, f'{taxon_id} is not a valid taxon ID.'
    print(f'Building protein tree for {id_to_names[taxon_id]} (ID: {taxon_id})...\n')
    run_protein_tree(user, password, taxon_id)
    print('Protein tree build done.')

if __name__ == '__main__':
  main()