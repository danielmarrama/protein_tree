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
    organism_ids = pd.DataFrame(result.fetchall(), columns=['organism_id']).astype(str)

    species = {} # map species taxon to all children organism IDs
    for organism_id in organism_ids['organism_id']:
      species_id = get_species_taxon_id(connection, organism_id)
      if species_id not in species:
        species[species_id] = []
      species[species_id].append(organism_id)
    
    species_list = [(species_id, '; '.join(organism_ids)) for species_id, organism_ids in species.items()]
    species_df = pd.DataFrame(species_list, columns=['species_taxon_id', 'organism_ids'])
    species_df.to_csv('species_test.csv', index=False)


def get_species_taxon_id(connection: Connection, organism_id: str) -> str:
  """
  Get the parent species taxon ID for an organism ID from the organism table.

  Args:
    connection: IEDB MySQL backend connection.
    organism_id: Organism ID.
  """
  path_query = """
              SELECT path
              FROM organism
              WHERE tax_id = :organism_id
              """
  result = connection.execute(text(path_query), {"organism_id": organism_id})
  path = result.fetchone()
  
  if path is None:
    return organism_id
  
  tax_ids = path[0].split(':')
  for tax_id in reversed(tax_ids):
    rank_query = """
                  SELECT rank
                  FROM organism
                  WHERE tax_id = :tax_id
                  """
    rank_result = connection.execute(text(rank_query), {"tax_id": tax_id})
    rank = rank_result.fetchone()[0]
    if rank == 'species':
      return tax_id

  return organism_id


if __name__ == '__main__':
  import argparse 

  parser = argparse.ArgumentParser(description='Update species data for protein tree.')
  parser.add_argument('-u', '--user', required=True, type=str, help='MySQL username.')
  parser.add_argument('-p', '--password', required=True, type=str, help='MySQL password.')

  args = parser.parse_args()
  update_species_data(args.user, args.password)