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
# * save source accessions from epitopes_df that are not in sources_df
# * use manual_parents.csv to override assigned genes
# * investigate a way to search all epitopes at once and make sure the assigned isoform is from the proper gene

def run_protein_tree(
  user: str, password: str, taxon_id: str, species_name: str, 
  species_df: pd.DataFrame, metrics_df: pd.DataFrame, all_taxa: str
) -> tuple:
  """
  Build protein tree for an IEDB species.
  
  Args:
    user: Username for IEDB MySQL connection.
    password: Password for IEDB MySQL connection.
    taxon_id: Taxon ID for the species to run protein tree.
    all_taxa: List of all children taxa for a species from the IEDB.
  """
  print('Getting epitopes and sources data...')

  Fetcher = DataFetcher(user, password)
  epitopes_df = Fetcher.get_epitopes(all_taxa)
  sources_df = Fetcher.get_sources(all_taxa)

  # if there are no epitopes or sources, return None
  if epitopes_df.empty or sources_df.empty:
    return None

  # or if the sources data is empty, return None
  if sources_df.dropna(subset=['Sequence']).empty:
    return None
  
  print('Done getting data.\n')

  # select best proteome for the species
  print('Getting the best proteome...')
  Selector = ProteomeSelector(taxon_id, species_name, species_df, metrics_df)
  print(f'Number of candidate proteomes: {Selector.num_of_proteomes}\n')

  proteome_data = Selector.select_best_proteome(epitopes_df)
  Selector.proteome_to_csv()
  
  print('Got the best proteome:')
  print(f'Proteome ID: {proteome_data[0]}')
  print(f'Proteome taxon: {proteome_data[1]}')
  print(f'Proteome type: {proteome_data[2]}\n')

  # assign genes to source antigens and parent proteins to epitopes
  print('Assigning genes to source antigens...')
  Assigner = GeneAssigner(taxon_id, species_name)
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
  idx = metrics_df['Species Taxon ID'] == int(taxon_id)
  
  metrics_df.loc[idx, 'Proteome ID'] = proteome_data[0]
  metrics_df.loc[idx, 'Proteome Taxon'] = proteome_data[1]
  metrics_df.loc[idx, 'Proteome Type'] = proteome_data[2]
  metrics_df.loc[idx, 'Source Count'] = assigner_data[0]
  metrics_df.loc[idx, 'Epitope Count'] = assigner_data[1]
  metrics_df.loc[idx, 'Successful Source Assignment'] = successful_source_assignment
  metrics_df.loc[idx, 'Successful Epitope Assignment'] = successful_epitope_assignment
  
  metrics_df.to_csv('metrics.csv', index=False)
  print('Done recording metrics.\n')

  return proteome_data, assigner_data


def build_tree():
  pass


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

  species_df = pd.read_csv('species.csv')
  metrics_df = pd.read_csv('metrics.csv')
  valid_taxon_ids = species_df['Species Taxon ID'].astype(str).tolist()

  # taxa and name mapppings
  all_taxa_map = dict(zip(species_df['Species Taxon ID'].astype(str), species_df['All Taxa']))
  species_id_to_name_map = dict(zip(species_df['Species Taxon ID'].astype(str), species_df['Species Name']))

  if all_species:
    for taxon_id in valid_taxon_ids:
      
      print(f'Building protein tree for {species_id_to_name_map[taxon_id]} (ID: {taxon_id})...\n')
      tree_data = run_protein_tree(
        user, password, taxon_id, species_id_to_name_map[taxon_id], 
        species_df, metrics_df, all_taxa_map[taxon_id]
      )

      if tree_data is None:
        print('No epitopes or sources found for this species. Skipping.')
        continue

      print('Protein tree build done.')

    print('All protein trees built.')

  else: # one species at a time
    assert taxon_id in valid_taxon_ids, f'{taxon_id} is not a valid taxon ID.'
    
    print(f'Building protein tree for {species_id_to_name_map[taxon_id]} (ID: {taxon_id})...\n')
    tree_data = run_protein_tree(
      user, password, taxon_id, species_id_to_name_map[taxon_id], 
      species_df, metrics_df, all_taxa_map[taxon_id]
    )    
    
    if tree_data == None:
      print('No epitopes or sources found for this species.')
      return
    
    print(f'Protein tree built for {species_id_to_name_map[taxon_id]} (ID: {taxon_id}).\n')

if __name__ == '__main__':
  main()