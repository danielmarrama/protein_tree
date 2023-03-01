#!/usr/bin/env python3

import argparse
import pandas as pd
from sqlalchemy import create_engine, text


def get_epitopes(taxon_id):
  """
  Get all epitopes for a species from the IEDB API.
  """
  pass
  # TODO: pull data from the MySQL database


def get_sources(taxon_id):
  """
  Get all source antigens for a species from the IEDB API.
  """
  pass


def main():
  # define command line args which will take in a taxon ID
  parser = argparse.ArgumentParser()
  
  parser.add_argument('-u', '--user', help='User for IEDB MySQL connection.')
  parser.add_argument('-p', '--password', help='User for IEDB MySQL connection.')
  
  args = parser.parse_args()
  user = args.user
  password = args.password

  sql_engine = create_engine(f'mysql://{user}:{password}@iedb-mysql.liai.org:33306/iedb_query')
  df = pd.DataFrame(sql_engine.connect().execute(text('SELECT * FROM source;')))
  print(df)

if __name__ == '__main__':
  main()
