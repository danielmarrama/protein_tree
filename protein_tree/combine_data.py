#!/usr/bin/env python3

import os
import pandas as pd
from pathlib import Path


def main():
  """Combine all source antigen and epitope assignments into two CSV files."""

  all_epitopes_df, all_sources_df = pd.DataFrame(), pd.DataFrame()
  
  data_path = Path(__file__).parent.parent / 'data'
  for species in (data_path / 'species').iterdir():
    epitopes_df = pd.read_csv(species / 'epitope_assginments.csv')
    sources_df = pd.read_csv(species / 'source_assignments.csv')
    
    # add taxon ID and species name columns
    taxon_id = species.name.split('-')[0]
    species_name = species.name.split('-')[1]

    epitopes_df['Taxon ID'] = taxon_id
    epitopes_df['Species Name'] = species_name

    sources_df['Taxon ID'] = taxon_id
    sources_df['Species Name'] = species_name

    # rearrange columns to put taxon ID and species name first
    epitopes_df = epitopes_df[['Taxon ID', 'Species Name'] + epitopes_df.columns[:-2].tolist()]
    sources_df = sources_df[['Taxon ID', 'Species Name'] + sources_df.columns[:-2].tolist()]

    all_epitopes_df = pd.concat([all_epitopes_df, epitopes_df])
    all_sources_df = pd.concat([all_sources_df, sources_df])

  all_epitopes_df.to_csv(data_path / 'all_epitope_assignments.csv', index=False)
  all_sources_df.to_csv(data_path / 'all_source_assignments.csv', index=False)

if __name__ == '__main__':
  main()