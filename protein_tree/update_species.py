#!/usr/bin/env python3

import os
import requests
import _pickle as pickle
import pandas as pd
from sqlalchemy import text
from sql_engine import create_sql_engine


def update_species_data(user, password):
  """Get all necessary species data and update species.csv."""
  sql_query = """
              SELECT object.organism2_id
              FROM epitope, object
              WHERE epitope.e_object_id = object.object_id
              AND object.object_sub_type IN ("Peptide from protein", "Discontinuous protein residues")
              UNION
              SELECT object.organism2_id
              FROM epitope, object
              WHERE epitope.related_object_id = object.object_id
              AND object.object_sub_type IN ("Peptide from protein", "Discontinuous protein residues")
              """
  # get all epitope data - map all IDs to species rank
  sql_engine = create_sql_engine(user, password)
  organism_df = pd.DataFrame(sql_engine.connect().execute(text(sql_query)))
  species_mapping = create_species_mapping(list(organism_df['organism2_id'].astype(str).unique()))

  # create dataframe with species ID, species name, and all taxa
  data = [(key, *value) for key, value in species_mapping.items()]
  df = pd.DataFrame(data, columns=['Lower Rank Taxon ID', 'Taxon Rank', 'Species Taxon ID', 'Species Name'])
  
  # group by Species Taxon ID and Species Name, aggregating the lower ranks
  grouped_df = df.groupby(['Species Taxon ID', 'Species Name'])['Lower Rank Taxon ID'].apply(lambda x: '; '.join(map(str, x))).reset_index()
  grouped_df = grouped_df.rename(columns={'Species Taxon ID': 'Taxon ID', 'Lower Rank Taxon ID': 'All Taxa'})
  grouped_df.to_csv(f'new_species_data.csv', index=False)

def create_species_mapping(taxon_ids):
  """Create mapping from lower rank taxonomy to species."""
  species_mapping = {}
  for taxon_id in taxon_ids:
    taxonomy_data = get_taxonomy_information(taxon_id)
    taxon_rank, species_id, species_name = get_species_data(taxon_id, taxonomy_data)
    species_mapping[taxon_id] = (taxon_rank, species_id, species_name)
  
  path = os.path.join(os.path.dirname(__file__), '../mappings/species_mapping.pickle')
  with open(path, 'wb') as f:
    pickle.dump(species_mapping, f)
  
  return species_mapping

def get_taxonomy_information(taxon_id):
  """Read in taxonomy information from UniProt."""
  url = f'https://rest.uniprot.org/taxonomy/{taxon_id}.json'
  response = requests.get(url)
  taxonomy_data = response.json()
  return taxonomy_data

def get_species_data(taxon_id, taxonomy_data):
  """Get the species rank data from the taxonomy information."""
  try:
    if taxonomy_data['rank'] in ['species', 'species group', 'genus', 'family']:
      return taxonomy_data['rank'], taxonomy_data['taxonId'], taxonomy_data['scientificName']
    else:
      for lineage_item in taxonomy_data['lineage']:
        if lineage_item['rank'] == 'species':
          return taxonomy_data['rank'], lineage_item['taxonId'], lineage_item['scientificName']
  except: # sometimes the taxon is not found
      return '', taxon_id, ''

if __name__ == '__main__':
  import argparse 

  parser = argparse.ArgumentParser(description='Update species data.')
  parser.add_argument('-u', '--user', required=True, type=str, help='MySQL username.')
  parser.add_argument('-p', '--password', required=True, type=str, help='MySQL password.')

  args = parser.parse_args()
  user = args.user
  password = args.password

  update_species_data(user, password)