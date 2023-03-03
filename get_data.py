#!/usr/bin/env python3

import argparse
import pandas as pd
from sqlalchemy import create_engine, text

# TODO: pull all data from the MySQL backend

class DataFetcher:
  """
  Fetch data from IEDB MySQL backend.
  """
  def __init__(self, user, password, taxon_id):
    self.sql_engine = create_engine(f'mysql://{user}:{password}@iedb-mysql.liai.org:33306/iedb_query')
    self.taxon_id = taxon_id

  def get_species(self):
    """
    Get all IEDB species.
    """
    return pd.read_csv('species.csv')

  def get_epitopes(self):
    """
    Get all epitopes for a species.
    """
    all_epitopes = pd.read_csv('snapshot_2022-12-20/upstream/epitopes.tsv', sep='\t')
    return all_epitopes[all_epitopes['Organism ID'].astype(str) == f'{taxon_id}.0']

  def get_sources(self):
    """
    Get all source antigens for a species.
    """
    return pd.DataFrame(self.sql_engine.connect().execute(text(f'SELECT * FROM source WHERE organism_id = {self.taxon_id};')))

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

  sql_engine = create_engine(f'mysql://{user}:{password}@iedb-mysql.liai.org:33306/iedb_query')

  sources_df = get_sources(sql_engine, taxon_id)

if __name__ == '__main__':
  main()
