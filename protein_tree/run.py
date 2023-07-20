#!/usr/bin/env python3

import os
import argparse
import pandas as pd
from pathlib import Path

from get_data import DataFetcher
from select_proteome import ProteomeSelector
from assign_gene_protein import GeneAndProteinAssigner


def run_protein_tree(
  taxon_id: int,
  all_taxa: list,
  species_name: str,
  is_vertebrate: bool,
  update_proteome: bool,
  epitopes_df: pd.DataFrame,
  sources_df: pd.DataFrame
) -> None:
  """Run all steps for the protein tree.
  
  Args:
    taxon_id: Taxon ID for the species to run protein tree.
    species_name: Name of the species to run protein tree.
    all_taxa: List of all children taxa for a species from the IEDB.
  """
  print(f'Building tree for {species_name} (ID: {taxon_id})...\n')
  
  # if there are no epitopes or sources, return None
  if epitopes_df.empty or sources_df.empty:
    print('No epitopes or sources found for this species.')
    return

  # update proteome if flag or if proteome doesn't exist
  proteome_file = Path(f'species/{taxon_id}-{species_name.replace(" ", "_")}/proteome.fasta')
  if update_proteome or not proteome_file.exists():
    print('Getting the best proteome...')
    Selector = ProteomeSelector(taxon_id, species_name)
    print(f'Number of candidate proteomes: {Selector.num_of_proteomes}\n')

    proteome_data = Selector.select_best_proteome(epitopes_df)
    Selector.proteome_to_csv()
    
    print('Got the best proteome:')
    print(f'Proteome ID: {proteome_data[0]}')
    print(f'Proteome taxon: {proteome_data[1]:.0f}')
    print(f'Proteome type: {proteome_data[2]}\n')

  # assign genes to source antigens and parent proteins to epitopes
  print('Assigning genes to source antigens...')
  Assigner = GeneAndProteinAssigner(taxon_id, species_name, is_vertebrate)
  assigner_data = Assigner.assign(sources_df, epitopes_df)
  print('Done assigning genes.\n')

  successful_source_assignment = (assigner_data[2] / assigner_data[0])*100
  successful_epitope_assignment = (assigner_data[3] / assigner_data[1])*100

  print(f'Number of sources: {assigner_data[0]}')
  print(f'Number of epitopes: {assigner_data[1]}')
  print(f'Successful source antigen assignments: {successful_source_assignment:.1f}%')
  print(f'Successful epitope assignments: {successful_epitope_assignment:.1f}%\n')

  # write data to metrics.csv
  metrics_df = pd.read_csv(Path(__file__).parent.parent / 'data' / 'metrics.csv')
  idx = metrics_df['Species Taxon ID'] == taxon_id
  
  if update_proteome:
    metrics_df.loc[idx, 'Proteome ID'] = proteome_data[0]
    metrics_df.loc[idx, 'Proteome Taxon'] = proteome_data[1]
    metrics_df.loc[idx, 'Proteome Type'] = proteome_data[2]
  
  metrics_df.loc[idx, 'Source Count'] = assigner_data[0]
  metrics_df.loc[idx, 'Epitope Count'] = assigner_data[1]
  metrics_df.loc[idx, 'Successful Source Assignment'] = successful_source_assignment
  metrics_df.loc[idx, 'Successful Epitope Assignment'] = successful_epitope_assignment
  
  metrics_df.to_csv('metrics.csv', index=False)


def build_tree_for_species(
  taxon_id: int,
  all_taxa_map: dict,
  species_name_map: dict,
  is_vertebrate_map: dict,
  update_proteome: bool,
  epitopes_df: pd.DataFrame,
  sources_df: pd.DataFrame
) -> None:
  """Build protein tree for a species.
  
  Args:
    taxon_id: Taxon ID for the species to run protein tree.
    all_taxa_map: Mapping of taxon ID to all children taxa.
    species_name_map: Mapping of taxon ID to species name.
    is_vertebrate_map: Mapping of taxon ID to is_vertebrate.
  """
  all_taxa = [int(taxon) for taxon in all_taxa_map[taxon_id].split(';')]
  species_name = species_name_map[taxon_id]
  is_vertebrate = is_vertebrate_map[taxon_id]

  run_protein_tree(
    taxon_id, all_taxa, species_name, is_vertebrate, 
    update_proteome, epitopes_df, sources_df
  )
  print(f'Protein tree built for {species_name} (ID: {taxon_id}).\n')


def main():
  parser = argparse.ArgumentParser()
  
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
  parser.add_argument(
    '-d', '--update_data',
    action='store_true',
    help='Pull the epitope and source tables from the IEDB backend.'
  )
  parser.add_argument(
    '-p', '--update_proteome',
    action='store_true',
    help='Update the proteome(s) to be used for the species.'
  )
  
  data_dir = Path(__file__).parent.parent / 'data'
  species_df = pd.read_csv(data_dir / 'species.csv')
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

  args = parser.parse_args()
  
  Fetcher = DataFetcher()
  if args.update_data:
    print('Getting all data...')
    Fetcher.get_all_data()
    print('All data written.')

  if args.all_species:
    for taxon_id in valid_taxon_ids[156:]:

      all_taxa = [int(taxon) for taxon in all_taxa_map[taxon_id].split(';')]
      epitopes_df = Fetcher.get_epitopes_for_species(all_taxa)
      sources_df = Fetcher.get_sources_for_species(all_taxa)
      build_tree_for_species(
        taxon_id, all_taxa_map, species_name_map, is_vertebrate_map, 
        args.update_proteome, epitopes_df, sources_df
      )

    print('All species complete.')

  else: # one species at a time
    taxon_id = args.taxon_id
    assert taxon_id in valid_taxon_ids, f'{taxon_id} is not a valid taxon ID.'

    all_taxa = [int(taxon) for taxon in all_taxa_map[taxon_id].split(';')]
    epitopes_df = Fetcher.get_epitopes_for_species(all_taxa)
    sources_df = Fetcher.get_sources_for_species(all_taxa)
    build_tree_for_species(
      taxon_id, all_taxa_map, species_name_map, is_vertebrate_map, 
      args.update_proteome, epitopes_df, sources_df
    )

if __name__ == '__main__':
  main()
