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
  df = pd.DataFrame(data, columns=['Lower Rank Taxon ID', 'Taxon Rank', 'Species Taxon ID', 'Species Name', 'Group', 'Vertebrate'])
  
  # group by Species, Taxon ID, and Species Name, aggregating the lower ranks
  grouped_df = df.groupby(['Species Taxon ID', 'Species Name', 'Group', 'Vertebrate'])['Lower Rank Taxon ID'].apply(lambda x: '; '.join(map(str, x))).reset_index()
  grouped_df = grouped_df.rename(columns={'Species Taxon ID': 'Taxon ID', 'Lower Rank Taxon ID': 'All Taxa'})
  grouped_df.to_csv(f'new_species_data.csv', index=False)

def create_species_mapping(taxon_ids):
  """Create mapping from lower rank taxonomy to species."""
  species_mapping = {}
  for taxon_id in taxon_ids:

    print(f'Getting data for: {taxon_id}...')

    species_data = get_species_data(taxon_id)

    if species_data is None:
      continue

    print(species_data)
    
    print(f'Taxon ID: {taxon_id}')
    print(f'Taxon Rank: {species_data["taxon_rank"]}')
    print(f'Species ID: {species_data["species_taxon_id"]}')
    print(f'Species Name: {species_data["species_name"]}')
    print(f'Group: {species_data["group"]}')
    print(f'Vertebrate: {species_data["vertebrate"]}')

    species_mapping[taxon_id] = (species_data['taxon_rank'], species_data['species_taxon_id'], species_data['species_name'], species_data['group'], species_data['vertebrate'])
  
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

def get_species_data(taxon_id):
  """Get the species rank data and group/vertebrate status from the taxonomy information."""
  taxonomy_data = get_taxonomy_information(taxon_id)

  if 'lineage' not in taxonomy_data:
    print(f"Taxon ID {taxon_id} not found in UniProt. Skipping.")
    return None
  
  group = ""
  is_vertebrate = 0
  species_data = {}

  # get group and vertebrate status
  for lineage_item in taxonomy_data["lineage"]:
    if lineage_item["rank"] == "superkingdom":
      group = lineage_item["scientificName"]
    if lineage_item["rank"] in ["class", "subclass", "no rank"]:
      if "Vertebrata" in lineage_item["scientificName"] or (
          "commonName" in lineage_item and "vertebrates" in lineage_item["commonName"]
      ):
        is_vertebrate = 1
        break

  # if the taxon is a species, species group, genus, or family, then the species data is the taxon data
  if taxonomy_data['rank'] in ['species', 'species group', 'genus', 'family']:
    species_data = {
      'taxon_rank': taxonomy_data['rank'],
      'species_taxon_id': taxonomy_data['taxonId'],
      'species_name': taxonomy_data['scientificName'],
      'group': group,
      'vertebrate': is_vertebrate,
    }
  else:
    for lineage_item in taxonomy_data['lineage']:
      if lineage_item['rank'] == 'species':
        species_data = {
            'taxon_rank': taxonomy_data['rank'],
            'species_taxon_id': lineage_item['taxonId'],
            'species_name': lineage_item['scientificName'],
            'group': group,
            'vertebrate': is_vertebrate,
        }
        break

  return species_data

if __name__ == '__main__':
  import argparse 

  parser = argparse.ArgumentParser(description='Update species data.')
  parser.add_argument('-u', '--user', required=True, type=str, help='MySQL username.')
  parser.add_argument('-p', '--password', required=True, type=str, help='MySQL password.')

  args = parser.parse_args()
  user = args.user
  password = args.password

  update_species_data(user, password)