#!/usr/bin/env python3

import pandas as pd
from sqlalchemy import create_engine, text


class DataFetcher:
  """
  Fetch data from IEDB MySQL backend.
  """
  def __init__(self, user, password, taxon_id, all_taxa):
    self.sql_engine = create_engine(f'mysql+mysqlconnector://{user}:{password}@iedb-mysql.liai.org:33306/iedb_query')
    self.taxon_id = taxon_id
    self.all_taxa = all_taxa.replace(';', ',')

  def get_epitopes(self):
    """
    Get all epitopes for a species.
    """
    sql_query = f'SELECT organism2_id, organism2_name, mol1_seq, '\
                f'mol2_name, mol2_accession '\
                f'FROM object '\
                f'WHERE object_sub_type = "Peptide from protein" '\
                f'AND organism2_id IN ({self.all_taxa});'
    columns = ['Organism ID', 'Organism_Name', 'Peptide', 'Source Name', 'Source Accession']
    return pd.DataFrame(self.sql_engine.connect().execute(text(sql_query)), columns=columns)

  def get_sources(self):
    """
    Get all source antigens for a species.
    """
    sql_query = f'SELECT accession, name, sequence, organism_id, organism_name '\
                f'FROM source WHERE organism_id IN ({self.all_taxa});'
    columns = ['Accession', 'Name', 'Sequence', 'Organism ID', 'Organism Name']
    return pd.DataFrame(self.sql_engine.connect().execute(text(sql_query)), columns=columns)

def main():
  import argparse
  import os

  parser = argparse.ArgumentParser()
  
  parser.add_argument('-u', '--user', required=True, help='User for IEDB MySQL connection.')
  parser.add_argument('-p', '--password', required=True, help='User for IEDB MySQL connection.')
  parser.add_argument('-t', '--taxon_id', required=True, help='Taxon ID for the species to pull data for.')
  
  args = parser.parse_args()
  user = args.user
  password = args.password
  taxon_id = args.taxon_id

  # read in IEDB species data
  species_df = pd.read_csv('species.csv')
  species_id_to_name_map = dict(zip(species_df['Taxon ID'].astype(str), species_df['Species Label']))
  all_taxa_map = dict(zip(species_df['Taxon ID'].astype(str), species_df['All Taxa']))

  # get epitopes and source antigens
  Fetcher = DataFetcher(user, password, taxon_id, all_taxa_map[taxon_id])
  epitopes_df = Fetcher.get_epitopes()
  sources_df = Fetcher.get_sources()

  # create directory for species and taxon ID
  species_path = f'species/{taxon_id}-{species_id_to_name_map[taxon_id].replace(" ", "_")}'
  os.makedirs(species_path, exist_ok=True)

  # write epitopes and source antigens to files
  epitopes_df.to_csv(f'{species_path}/epitopes.csv', index=False)
  sources_df.to_csv(f'{species_path}/sources.csv', index=False)

if __name__ == '__main__':
  main()
