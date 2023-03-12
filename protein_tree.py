#!/usr/bin/env python3

import os
import argparse
import pandas as pd
import numpy as np

from get_data import DataFetcher
from select_proteome import ProteomeSelector
from assign_genes import GeneAssigner

# TODO:
# - create UniProt ID to gene symbol map using proteome.csv
# - create a .txt file with the tree structure of gene --> relevant isoforms
# - either update PEPMatch to search discontinous epitopes or write a new function to do it

def run_protein_tree(user, password, all_species, taxon_id, species_name, all_taxa):
  """
  Build protein tree for an IEDB species.
  """

  # get epitopes and sources data from MySQL backend
  print('Getting epitopes and sources data...')

  Fetcher = DataFetcher(user, password, taxon_id, all_taxa)
  epitopes_df = Fetcher.get_epitopes()
  sources_df = Fetcher.get_sources()

  # if there are no epitopes or sources, return the proteome data that exists in metrics.csv and then zeros for the other tuple
  if epitopes_df.empty or sources_df.empty:
    return None
  
  print('Done getting data.\n')

  # select best proteome for the species
  print('Getting the best proteome...')
  Selector = ProteomeSelector(taxon_id)
  print(f'Number of candidate proteomes: {Selector.num_of_proteomes}\n')

  proteome_data = Selector.select_proteome(epitopes_df)
  Selector.proteome_to_csv()
  
  print('Got the best proteome:')
  print(f'Proteome ID: {proteome_data[0]}')
  print(f'Proteome taxon: {proteome_data[1]}')
  print(f'Proteome type: {proteome_data[2]}\n')

  # assign genes to source antigens and parent proteins to epitopes
  print('Assigning genes to source antigens...')
  Assigner = GeneAssigner(taxon_id)
  assigner_data = Assigner.assign_genes(sources_df, epitopes_df)
  print('Done assigning genes.\n')

  print(f'Number of sources: {assigner_data[0]}')
  print(f'Number of epitopes: {assigner_data[4]}')
  print(f'Number of sources missing sequences: {assigner_data[1]}')
  print(f'Number of sources with no BLAST matches: {assigner_data[2]}')
  print(f'Number of sources with BLAST matches: {assigner_data[3]}')
  print(f'Number of epitopes with a match: {assigner_data[5]}')
  print(f'Successful gene assignments: {(assigner_data[3] / assigner_data[0])*100:.1f}%')
  print(f'Successful parent assignments: {(assigner_data[5] / assigner_data[4])*100:.1f}%\n')

  # write data to metrics.csv
  print('Recording metrics...') 
  metrics_df = pd.read_csv('metrics.csv')
  metrics_df.loc[metrics_df['Taxon ID'] == int(taxon_id), 'Proteome ID'] = proteome_data[0]
  metrics_df.loc[metrics_df['Taxon ID'] == int(taxon_id), 'Proteome Taxon'] = proteome_data[1]
  metrics_df.loc[metrics_df['Taxon ID'] == int(taxon_id), 'Proteome Type'] = proteome_data[2]
  metrics_df.loc[metrics_df['Taxon ID'] == int(taxon_id), '# of Sources'] = assigner_data[0]
  metrics_df.loc[metrics_df['Taxon ID'] == int(taxon_id), '# of Epitopes'] = assigner_data[4]
  metrics_df.loc[metrics_df['Taxon ID'] == int(taxon_id), '# of Sources Missing Sequences'] = assigner_data[1]
  metrics_df.loc[metrics_df['Taxon ID'] == int(taxon_id), '# of Sources with No BLAST Matches'] = assigner_data[2]
  metrics_df.loc[metrics_df['Taxon ID'] == int(taxon_id), '# of Sources with BLAST Matches'] = assigner_data[3]
  metrics_df.loc[metrics_df['Taxon ID'] == int(taxon_id), '# of Epitopes with a Match'] = assigner_data[5]
  
  # calculate percentages of gene and parent assignment success
  metrics_df.loc[metrics_df['Taxon ID'] == int(taxon_id), '% Successful Gene Assignments'] = (assigner_data[3] / assigner_data[0])*100
  metrics_df.loc[metrics_df['Taxon ID'] == int(taxon_id), '% Successful Parent Assignments'] = (assigner_data[5] / assigner_data[4])*100

  metrics_df.to_csv('metrics.csv', index=False)
  print('Done recording metrics.\n')

  return (proteome_data, assigner_data)

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

  run_protein_tree(user, password, all_species, taxon_id)

  # read in IEDB species data and read in metrics data to write into
  species_df = pd.read_csv('species.csv')
  metrics_df = pd.read_csv('metrics.csv')
  valid_taxon_ids = species_df['Taxon ID'].astype(str).tolist()

  # dicts for mapping taxon IDs to all their taxa and their names
  all_taxa_map = dict(zip(species_df['Taxon ID'].astype(str), species_df['All Taxa']))
  species_id_to_name_map = dict(zip(species_df['Taxon ID'].astype(str), species_df['Species Label']))

  # run protein tree for all IEDB species 
  if all_species:
    for t_id in valid_taxon_ids:
      print(f'Building protein tree for {species_id_to_name_map[t_id]} (ID: {t_id})...\n')
      tree_data = run_protein_tree(user, password, t_id, species_id_to_name_map[t_id], all_taxa_map[t_id])

      if tree_data is None:
        print('No epitopes or sources found for this species. Skipping.')
        continue

      print('Protein tree build done.')

    print('All protein trees built.')

  # or one species at a time
  else:
    assert taxon_id in valid_taxon_ids, f'{taxon_id} is not a valid taxon ID.'
    print(f'Building protein tree for {species_id_to_name_map[taxon_id]} (ID: {taxon_id})...\n')
    tree_data = run_protein_tree(user, password, taxon_id, species_id_to_name_map[taxon_id], all_taxa_map[taxon_id])
    
    if tree_data == None:
      print('No epitopes or sources found for this species.')
      return
    
    print(f'Protein tree built for {species_id_to_name_map[taxon_id]} (ID: {taxon_id}).\n')

if __name__ == '__main__':
  main()