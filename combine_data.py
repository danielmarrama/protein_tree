#!/usr/bin/env python3

import os
import pandas as pd

# loop through every species dir in species/ and combine all epitopes.csv files into one file
# and all sources.csv files into one file

def main():
  all_epitopes_df, all_sources_df = pd.DataFrame(), pd.DataFrame()
  for species in os.listdir('species/'):
    try:
      epitopes_df = pd.read_csv(f'species/{species}/epitopes.csv')
      sources_df = pd.read_csv(f'species/{species}/sources.csv')
    except FileNotFoundError:
      continue

    # add taxon ID and species name columns
    epitopes_df['Taxon ID'] = species.split('-')[0]
    epitopes_df['Species Name'] = species.split('-')[1]

    sources_df['Taxon ID'] = species.split('-')[0]
    sources_df['Species Name'] = species.split('-')[1]

    # rearrange columns to put taxon ID and species name first
    epitopes_df = epitopes_df[['Taxon ID', 'Species Name'] + epitopes_df.columns[:-2].tolist()]
    sources_df = sources_df[['Taxon ID', 'Species Name'] + sources_df.columns[:-2].tolist()]

    all_epitopes_df = pd.concat([all_epitopes_df, epitopes_df])
    all_sources_df = pd.concat([all_sources_df, sources_df])

  all_epitopes_df.to_csv('all_epitopes.csv', index=False)
  all_sources_df.to_csv('all_sources.csv', index=False)

if __name__ == '__main__':
  main()