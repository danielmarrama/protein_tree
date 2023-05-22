#!/usr/bin/env python3

import requests
import _pickle as pickle
import pandas as pd
from pathlib import Path
from sqlalchemy import text
from sqlalchemy.engine import Connection
from protein_tree.sql_engine import create_sql_engine


def update_species_data(user: str, password: str) -> None:
  """
  Get all organism IDs for all epitope data we need for protein tree. Then,
  get the species taxon ID for each organism and update the species data file.

  Args:
    user: IEDB MySQL backend username.
    password: IEDB MySQL backend password.
  """
  sql_query = """
              SELECT object.organism2_id, object.organism2_name
              FROM epitope, object
              WHERE epitope.e_object_id = object.object_id
              AND object.object_sub_type IN ("Peptide from protein", "Discontinuous protein residues")
              UNION
              SELECT object.organism2_id, object.organism2_name
              FROM epitope, object
              WHERE epitope.related_object_id = object.object_id
              AND object.object_sub_type IN ("Peptide from protein", "Discontinuous protein residues")
              """
  sql_engine = create_sql_engine(user, password)
  with sql_engine.connect() as connection:
    
    result = connection.execute(text(sql_query))
    organisms = pd.DataFrame(result.fetchall())
    organisms.columns = ['organism_id', 'organism_name']
    organisms = organisms[organisms['organism_id'].notna()] # remove null organism IDs

    species = {} # map species taxon to all children organism IDs, is_vertebrate, and group
    for organism in organisms.itertuples():

      species_data = get_species_data(connection, str(organism.organism_id), organism.organism_name)

      if species_data[0] not in species: # add new species
        species[species_data[0]] = {
          'species_name': species_data[1], 
          'all_taxa': [], 
          'is_vertebrate': species_data[2], 
          'group': species_data[3]
        }
      species[species_data[0]]['all_taxa'].append(str(organism.organism_id))

    species_df = pd.DataFrame.from_dict(species, orient='index')
    species_df.index.name = 'species_taxon_id'
    species_df.reset_index(inplace=True)
    species_df['all_taxa'] = species_df['all_taxa'].apply(lambda x: ';'.join(x))
    species_df.to_csv('species.csv', index=False)


def get_species_data(connection: Connection, organism_id: str, organism_name: str) -> tuple:
  """
  Using the organism taxon ID from the IEDB backend, get the parent species
  taxon ID, the superkingdom, and whether or not the species is a vertebrate. 

  Args:
    connection: IEDB MySQL backend connection.
    organism_id: Organism taxon ID.
    organism_name: Organism name.
  """
  path_query = """
              SELECT path
              FROM organism
              WHERE tax_id = :organism_id
              """
  result = connection.execute(text(path_query), {"organism_id": organism_id})
  path = result.fetchone()
  
  if path is None:
    return (organism_id, organism_name, False, 'Other')
  
  tax_ids = path[0].split(':')

  # assign is_vertebrate
  is_vertebrate = True if '7742' in tax_ids else False # Vertebrata taxon ID
  
  # assign superkingdom
  groups = {'2': 'Bacteria', '2157': 'Archaea', '2759': 'Eukaryota', '10239': 'Viruses'} 
  group = 'Other'
  for key in groups.keys(): 
    if key in tax_ids:
      group = groups[key]
      break

  # assign species taxon ID and species name
  species_id = organism_id
  species_name = organism_name
  for tax_id in reversed(tax_ids):
    rank_query = """
                  SELECT rank, organism_name
                  FROM organism
                  WHERE tax_id = :tax_id
                  """
    rank_result = connection.execute(text(rank_query), {"tax_id": tax_id})
    rank, name = rank_result.fetchone()

    if rank == 'species':
      species_id = tax_id
      species_name = name
      break
  
  return (species_id, species_name, is_vertebrate, group)


if __name__ == '__main__':
  import argparse 

  parser = argparse.ArgumentParser(description='Update species data for protein tree.')
  parser.add_argument('-u', '--user', required=True, type=str, help='MySQL username.')
  parser.add_argument('-p', '--password', required=True, type=str, help='MySQL password.')

  args = parser.parse_args()
  update_species_data(args.user, args.password)