#!/usr/bin/env python3

import os
import argparse
import pandas as pd
import numpy as np

from get_data import DataFetcher
from select_proteome import ProteomeSelector
from assign_genes import GeneAssigner

# MAIN TODO:
# * add step to handle allergens differntly by assigning label from IUIS
# * run ARC before BLAST
# * add step to override gene and parent assignments with manual assignments

# smaller TODO:
# * use manual_parents.csv to override assigned genes
# * investigate a way to search all epitopes at once and make sure the assigned isoform is from the proper gene

def run_protein_tree(
  user: str,
  password: str,
  taxon_id: int,
  species_name: str,
  all_taxa: str,
  is_vertebrate: bool
) -> int:
  """Run all steps for the protein tree.
  
  Args:
    user: Username for IEDB MySQL connection.
    password: Password for IEDB MySQL connection.
    taxon_id: Taxon ID for the species to run protein tree.
    species_name: Name of the species to run protein tree.
    all_taxa: List of all children taxa for a species from the IEDB.
  """
  species_df = pd.read_csv('species.csv') # all species and their taxon IDs

  print('Getting epitopes and sources data...')
  Fetcher = DataFetcher(user, password)
  epitopes_df = Fetcher.get_epitopes(all_taxa)
  sources_df = Fetcher.get_sources(all_taxa)
  print('Done getting data.\n')
  
  # if there are no epitopes or sources, return None
  if epitopes_df.empty or sources_df.empty:
    return None

  print('Getting the best proteome...')
  Selector = ProteomeSelector(taxon_id, species_name)
  print(f'Number of candidate proteomes: {Selector.num_of_proteomes}\n')

  proteome_data = Selector.select_best_proteome(epitopes_df)
  Selector.proteome_to_csv()
  
  print('Got the best proteome:')
  print(f'Proteome ID: {proteome_data[0]}')
  print(f'Proteome taxon: {proteome_data[1]}')
  print(f'Proteome type: {proteome_data[2]}\n')

  # assign genes to source antigens and parent proteins to epitopes
  print('Assigning genes to source antigens...')
  Assigner = GeneAssigner(taxon_id, species_name, is_vertebrate)
  assigner_data = Assigner.assign_genes(sources_df, epitopes_df)
  print('Done assigning genes.\n')

  successful_source_assignment = (assigner_data[2] / assigner_data[0])*100
  successful_epitope_assignment = (assigner_data[3] / assigner_data[1])*100

  print(f'Number of sources: {assigner_data[0]}')
  print(f'Number of epitopes: {assigner_data[1]}')
  print(f'Number of sources with a match: {assigner_data[2]}')
  print(f'Number of epitopes with a match: {assigner_data[3]}')
  print(f'Successful source antigen assignments: {successful_source_assignment:.1f}%')
  print(f'Successful epitope assignments: {successful_epitope_assignment:.1f}%\n')

  # write data to metrics.csv
  print('Recording metrics...') 
  metrics_df = pd.read_csv('metrics.csv')
  idx = metrics_df['Species Taxon ID'] == taxon_id
  
  metrics_df.loc[idx, 'Proteome ID'] = proteome_data[0]
  metrics_df.loc[idx, 'Proteome Taxon'] = proteome_data[1]
  metrics_df.loc[idx, 'Proteome Type'] = proteome_data[2]
  metrics_df.loc[idx, 'Source Count'] = assigner_data[0]
  metrics_df.loc[idx, 'Epitope Count'] = assigner_data[1]
  metrics_df.loc[idx, 'Successful Source Assignment'] = successful_source_assignment
  metrics_df.loc[idx, 'Successful Epitope Assignment'] = successful_epitope_assignment
  
  metrics_df.to_csv('metrics.csv', index=False)
  print('Done recording metrics.\n')

  return 1


def build_tree_for_species(
  user: str,
  password: str,
  taxon_id: int,
  all_taxa_map: dict,
  species_name_map: dict,
  is_vertebrate_map: dict
) -> None:
  """Build protein tree for a species.
  
  Args:
    user: Username for IEDB MySQL connection.
    password: Password for IEDB MySQL connection.
    taxon_id: Taxon ID for the species to run protein tree.
    all_taxa_map: Mapping of taxon ID to all children taxa.
    species_name_map: Mapping of taxon ID to species name.
    is_vertebrate_map: Mapping of taxon ID to is_vertebrate.
  """
  all_taxa = all_taxa_map[taxon_id]
  species_name = species_name_map[taxon_id]
  is_vertebrate = is_vertebrate_map[taxon_id]

  run_protein_tree(
    user, password, taxon_id, species_name, all_taxa, is_vertebrate
  )
  print(f'Protein tree built for {species_name} (ID: {taxon_id}).\n')


def main():
  parser = argparse.ArgumentParser()
  
  parser.add_argument(
    '-u', '--user',
    required=True,
    help='User for IEDB MySQL connection.'
  )
  parser.add_argument(
    '-p', '--password',
    required=True,
    help='Password for IEDB MySQL connection.'
  )
  parser.add_argument(
    '-a', '--all_species', 
    action='store_true',
    help='Build protein tree for all IEDB species.'
  )
  parser.add_argument(
    '-t', '--taxon_id', 
    type=int, 
    help='Taxon ID for the species to run protein tree.'
  )
  
  args = parser.parse_args()
  user = args.user
  password = args.password

  species_df = pd.read_csv('species.csv') # all species and their taxon IDs
  valid_taxon_ids = species_df['Species Taxon ID'].tolist()

  # taxa, species name, and is_vertebrate mapppings
  all_taxa_map = dict(
    zip(
      species_df['Species Taxon ID'],
      species_df['All Taxa']
    )
  )
  species_name_map = dict(
    zip(
      species_df['Species Taxon ID'],
      species_df['Species Name']
    )
  )
  is_vertebrate_map = dict(
    zip(
      species_df['Species Taxon ID'],
      species_df['Is Vertebrate']
    )
  )

  if args.all_species:
    for taxon_id in valid_taxon_ids:
      build_tree_for_species(
        user, password, taxon_id, all_taxa_map, species_name_map, is_vertebrate_map
      )
    print('All protein trees built.')

  else: # one species at a time
    taxon_id = args.taxon_id
    assert taxon_id in valid_taxon_ids, f'{taxon_id} is not a valid taxon ID.'
    build_tree_for_species(
      user, password, taxon_id, all_taxa_map, species_name_map, is_vertebrate_map
    )

if __name__ == '__main__':
  main()