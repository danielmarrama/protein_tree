#!/usr/bin/env python3

import os
import argparse
import pandas as pd

from get_data import DataFetcher
from select_proteome import ProteomeSelector
from assign_genes import GeneAssigner

# TODO:
# - when running all species, update species.csv with proteome ID, proteome taxon, and proteome type
#

def run_protein_tree(user, password, taxon_id, species_name, all_taxa):
  """
  Build protein tree for a species.
  """
  # make species folder if it doesn't exist
  os.makedirs(f'species/{taxon_id}-{species_name.replace(" ", "_")}', exist_ok=True)

  print('Getting epitopes and sources data...')

  Fetcher = DataFetcher(user, password, taxon_id, all_taxa)
  epitopes_df = Fetcher.get_epitopes()
  sources_df = Fetcher.get_sources()
  
  print('Done getting data.\n')

  # print('Getting the best proteome...')
  # Selector = ProteomeSelector(taxon_id)
  # proteome_id, proteome_taxon, proteome_type = Selector.select_proteome(epitopes_df)
  # Selector.proteome_to_csv()
  
  # print(f'Number of candidate proteomes: {Selector.num_of_proteomes}\n')
  # print('Got the best proteome:')
  # print(f'Proteome ID: {proteome_id}')
  # print(f'Proteome taxon: {proteome_taxon}')
  # print(f'Proteome type: {proteome_type}')

  print('\nAssigning genes to source antigens...\n')
  Assigner = GeneAssigner(taxon_id)
  Assigner.assign_genes(sources_df, epitopes_df)
  
  # Assigner.assign_parents()

def main():
  parser = argparse.ArgumentParser()
  
  parser.add_argument('-u', '--user', required=True, help='User for IEDB MySQL connection.')
  parser.add_argument('-p', '--password', required=True, help='Password for IEDB MySQL connection.')
  parser.add_argument('-a', '--all_species', action='store_true', help='Build protein tree for all IEDB species.')
  parser.add_argument('-t', '--taxon_id', help='Taxon ID for the species to run protein tree.')
  
  args = parser.parse_args()
  user = args.user
  password = args.password
  all_species = args.all_species
  taxon_id = args.taxon_id

  # read in IEDB species data
  species_df = pd.read_csv('species.csv')
  valid_taxon_ids = species_df['Taxon ID'].astype(str).tolist()

  # dicts for mapping taxon IDs to all their taxa and their names
  all_taxa_map = dict(zip(species_df['Taxon ID'].astype(str), species_df['All Taxa']))
  species_id_to_name_map = dict(zip(species_df['Taxon ID'].astype(str), species_df['Species Label']))

  # make species folder if it doesn't exist
  os.makedirs('species', exist_ok=True)

  # run protein tree for all IEDB species 
  if all_species:
    for t_id in valid_taxon_ids:
      print(f'Building protein tree for {species_id_to_name_map[t_id]} (ID: {t_id})...\n')
      run_protein_tree(user, password, t_id, species_id_to_name_map[t_id], all_taxa_map[t_id])
      print('Protein tree build done.')

  # or one species at a time
  else:
    assert taxon_id in valid_taxon_ids, f'{taxon_id} is not a valid taxon ID.'
    print(f'Building protein tree for {species_id_to_name_map[taxon_id]} (ID: {taxon_id})...\n')
    run_protein_tree(user, password, taxon_id, species_id_to_name_map[taxon_id], all_taxa_map[taxon_id])
    print('Protein tree build done.')

if __name__ == '__main__':
  main()