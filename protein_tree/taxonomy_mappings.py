#!/usr/bin/env python3

import os
import requests
import _pickle as pickle


def get_taxonomy_information(taxon_id):
  url = f'https://rest.uniprot.org/taxonomy/{taxon_id}.json'
  response = requests.get(url)
  taxonomy_data = response.json()
  return taxonomy_data

def get_species_data(taxon_id, taxonomy_data):
  try:
    if taxonomy_data['rank'] in ['species', 'species group', 'genus', 'family']:
      return taxonomy_data['rank'], taxonomy_data['taxonId'], taxonomy_data['scientificName']
    else:
      for lineage_item in taxonomy_data['lineage']:
        if lineage_item['rank'] == 'species':
          return taxonomy_data['rank'], lineage_item['taxonId'], lineage_item['scientificName']
  except:
      return '', taxon_id, ''

def get_species_children(taxon_id):
  url = f'https://rest.uniprot.org/taxonomy/stream?format=json&query=(ancestor:{taxon_id})'
  response = requests.get(url)
  taxonomy_data = response.json()
  return taxonomy_data

def create_species_mapping(taxon_ids):
  species_mapping = {}
  for taxon_id in taxon_ids:
    taxonomy_data = get_taxonomy_information(taxon_id)
    taxon_rank, species_id, species_name = get_species_data(taxon_id, taxonomy_data)
    species_mapping[taxon_id] = (taxon_rank, species_id, species_name)

  path = os.path.join(os.path.dirname(__file__), '../mappings/species_mapping.pickle')
  with open(path, 'wb') as f:
    pickle.dump(species_mapping, f)

  return species_mapping