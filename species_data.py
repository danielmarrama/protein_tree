#!/usr/bin/env python3

import os
import pandas as pd

if __name__ == '__main__':
  df = pd.DataFrame()
  for path, _, files in os.walk('build/species'):
    for name in files:
      if name == 'species-data.tsv':
        species_data_df = pd.read_csv(os.path.join(path, name), sep='\t')
        species_data_df['Species Taxon ID'] = path.split('/')[-1]
        df = pd.concat([df, species_data_df])

  df.reset_index(inplace=True, drop=True)
  df.insert(0, 'Species Taxon ID', df.pop('Species Taxon ID')) # move taxon ID to first
  df.to_csv('build/species-data.tsv', sep='\t', index=False)