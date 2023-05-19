#!/usr/bin/env python3

import pandas as pd
from pathlib import Path
from sqlalchemy import text
from sql_engine import create_sql_engine


class DataFetcher:
  def __init__(self, user: str, password: str) -> None:
    self.sql_engine = create_sql_engine(user, password) # private so there's no exposure to the backend

  def get_epitopes(self, all_taxa: list) -> pd.DataFrame:
    """
    Get all epitopes for a species including all children taxa.

    Args:
      all_taxa: List of all children taxa for a species from the IEDB.
    """
    all_taxa = all_taxa.replace(';', ',')
    sql_query1 = f"""
                  SELECT object.mol1_seq, object.mol2_name, object.mol2_accession
                  FROM epitope, object
                  WHERE epitope.e_object_id = object.object_id
                  AND object.object_sub_type IN ("Peptide from protein", "Discontinuous protein residues")
                  AND object.organism2_id IN ({all_taxa});
                  """
    sql_query2 = f"""
                  SELECT object.mol1_seq, object.mol2_name, object.mol2_accession
                  FROM epitope, object
                  WHERE epitope.related_object_id = object.object_id
                  AND object.object_sub_type IN ("Peptide from protein", "Discontinuous protein residues")
                  AND object.organism2_id IN ({all_taxa});
                  """
    
    with self.sql_engine.connect() as connection:
      print('Fetching epitopes...')
      result1 = connection.execute(text(sql_query1))
      result2 = connection.execute(text(sql_query2))

      df1 = pd.DataFrame(result1.fetchall(), columns=['Sequence', 'Source Name', 'Source Accession'])
      df2 = pd.DataFrame(result2.fetchall(), columns=['Sequence', 'Source Name', 'Source Accession'])
      epitopes_df = pd.concat([df1, df2], ignore_index=True)
    
    return epitopes_df

  def get_sources(self, all_taxa: list) -> pd.DataFrame:
    """
    Get all source antigens for a species including all children taxa.

    Args:
      all_taxa: List of all children taxa for a species from the IEDB.
    """
    all_taxa = all_taxa.replace(';', ',')
    sql_query = f'SELECT source_id, accession, name, sequence '\
                f'FROM source WHERE organism_id IN ({all_taxa});'

    with self.sql_engine.connect() as connection:
      print('Fetching source antigens...')
      result = connection.execute(text(sql_query))
      sources_df = pd.DataFrame(result.fetchall(), columns=['Source ID', 'Accession', 'Name', 'Sequence'])
    
    return sources_df

def main():
  import argparse
  import os

  parser = argparse.ArgumentParser()
  
  parser.add_argument('-u', '--user', required=True, help='User for IEDB MySQL connection.')
  parser.add_argument('-p', '--password', required=True, help='User for IEDB MySQL connection.')
  parser.add_argument('-a', '--all_species', action='store_true', help='Build protein tree for all IEDB species.')
  parser.add_argument('-t', '--taxon_id', type=int, help='Taxon ID for species to pull data for.')
  
  args = parser.parse_args()

  Fetcher = DataFetcher(args.user, args.password)
  if args.all_species:
    Fetcher.get_species_data()
  else:
    taxon_id = args.taxon_id

    parent_dir = Path(__file__).parent
    species_df = pd.read_csv(parent_dir / '..' / 'species.csv')
    
    # create maps for taxon ID to species name and all taxa
    species_name_map = dict(zip(species_df['species_taxon_id'], species_df['species_name']))
    all_taxa_map = dict(zip(species_df['species_taxon_id'], species_df['all_taxa']))

    # get epitopes and source antigens
    epitopes_df = Fetcher.get_epitopes(all_taxa_map[taxon_id])
    sources_df = Fetcher.get_sources(all_taxa_map[taxon_id])

    # create directory for species and taxon ID
    species_path = parent_dir / '..' / 'species' / f'{taxon_id}-{species_name_map[taxon_id].replace(" ", "_")}'
    os.makedirs(species_path, exist_ok=True)

    # write epitopes and source antigens to files
    epitopes_df.to_csv(species_path / 'epitopes.csv', index=False)
    sources_df.to_csv(species_path / 'sources.csv', index=False)

if __name__ == '__main__':
  main()
