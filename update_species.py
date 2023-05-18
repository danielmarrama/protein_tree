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
  sql_engine = create_sql_engine(user, password)
  with sql_engine.connect() as connection:
    result = connection.execute(text(sql_query))
    organism_ids = pd.DataFrame(result.fetchall(), columns=['organism_id'])
    organism_ids = organism_ids[organism_ids['organism_id'].notna()] # remove null organism IDs

    species = {} # map species taxon to all children organism IDs, is_vertebrate, and group
    for organism_id in organism_ids['organism_id']:
      print(organism_id)
      species_data = get_species_data(connection, organism_id)
      
      if species_data[0] not in species:
        species[species_data[0]] = {'organism_ids': [], 'is_vertebrate': species_data[1], 'group': species_data[2]}
      
      species[species_data[0]]['organism_ids'].append(organism_id)
    
    species_list = [(species_id, '; '.join(data["organism_ids"]), data["is_vertebrate"], data["group"]) 
                    for species_id, data in species.items()] # convert dict to tuples for dataframe
    
    species_df = pd.DataFrame(species_list, columns=['Species Taxon ID', 'All Taxa', 'Is Vertebrate', 'Group'])
    species_df.to_csv('species.csv', index=False)


def get_species_data(connection: Connection, organism_id: str) -> tuple:
  """
  Using the organism taxon ID from the IEDB backend, get the parent species
  taxon ID, the superkingdom, and whether or not the species is a vertebrate. 

  Args:
    connection: IEDB MySQL backend connection.
    organism_id: Organism taxon ID.
  """
  path_query = """
              SELECT path
              FROM organism
              WHERE tax_id = :organism_id
              """
  result = connection.execute(text(path_query), {"organism_id": str(organism_id)})
  path = result.fetchone()
  
  if path is None:
    return (organism_id, False, 'Other')
  
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

  # assign species taxon ID
  species_id = organism_id
  for tax_id in reversed(tax_ids):
    rank_query = """
                  SELECT rank
                  FROM organism
                  WHERE tax_id = :tax_id
                  """
    rank_result = connection.execute(text(rank_query), {"tax_id": tax_id})
    rank = rank_result.fetchone()[0]

    if rank == 'species':
      species_id = tax_id
      break

  return (species_id, is_vertebrate, group)


if __name__ == '__main__':
  import argparse 

  parser = argparse.ArgumentParser(description='Update species data for protein tree.')
  parser.add_argument('-u', '--user', required=True, type=str, help='MySQL username.')
  parser.add_argument('-p', '--password', required=True, type=str, help='MySQL password.')

  args = parser.parse_args()
  update_species_data(args.user, args.password)