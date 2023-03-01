#!/usr/bin/env python3

import requests

def get_epitopes(taxon_id):
  """
  Get all epitopes for a species from the IEDB API.
  """

  # TODO: pull data from the MySQL database

  url = f'https://query-api.iedb.org/epitope_search'

def get_sources(taxon_id):
  """
  Get all source antigens for a species from the IEDB API.
  """
  url = f'https://query-api.iedb.org/parent_protein'


def main():
  # define command line args which will take in a taxon ID
  parser = argparse.ArgumentParser(description='Select the taxon ID to get the epitopes and source antigens.')
  parser.add_argument('taxon_id', help='Taxon ID for the species.')
  args = parser.parse_args()
  taxon_id = args.taxon_id

  # get all epitopes for the species
  get_epitopes(taxon_id)

  # get all source antigens for the species
  get_sources(taxon_id)

if __name__ == '__main__':
  main()
