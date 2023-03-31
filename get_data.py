#!/usr/bin/env python3

import pandas as pd
from sql_engine import create_sql_engine

class DataFetcher:
  """
  Fetch data from IEDB MySQL backend.
  """
  def __init__(self, user, password):
    self.sql_engine = create_sql_engine(user, password) # private so there's no exposure to the backend

  def get_species(self):
    sql_query = """
                SELECT object.*
                FROM iedb_curation.epitope epitope, iedb_curation.object object
                WHERE epitope.epitope_id = object.object_id
                AND object.object_sub_type IN ("Peptide from protein", "Discontinuous protein residues")
                UNION
                SELECT object.*
                FROM iedb_curation.epitope epitope, iedb_curation.object object
                WHERE epitope.related_object_id = object.object_id
                AND object.object_sub_type IN ("Peptide from protein", "Discontinuous protein residues")
                """
    # read in all epitope data from IEDB
    df = pd.DataFrame(self.sql_engine.connect().execute(text(sql_query)))
    
    # zip species IDs and names into a dictionary and save to .csv file using dataframe
    pd.DataFrame.from_dict(dict(zip(df['organism2_id'].astype(int), df['organism2_name'])), 'index').reset_index().to_csv('species2.csv')

  def get_epitopes(self, all_taxa):
    """
    Get all epitopes for a species.
    """
    all_taxa = all_taxa.replace(';', ',')
    sql_query = f'SELECT mol1_seq, mol2_name, mol2_accession '\
                f'FROM object '\
                f'WHERE object_sub_type = "Peptide from protein" '\
                f'AND organism2_id IN ({self.all_taxa});'
    columns = ['Sequence', 'Source Name', 'Source Accession']
    epitopes_df = pd.DataFrame(self.sql_engine.connect().execute(text(sql_query)), columns=columns)
    return epitopes_df.drop_duplicates(subset=['Sequence', 'Source Name', 'Source Accession'])

  def get_sources(self, all_taxa):
    """
    Get all source antigens for a species.
    """
    all_taxa = all_taxa.replace(';', ',')
    sql_query = f'SELECT source_id, accession, name, sequence '\
                f'FROM source WHERE organism_id IN ({self.all_taxa});'
    columns = ['Source ID', 'Accession', 'Name', 'Sequence']
    sources_df = pd.DataFrame(self.sql_engine.connect().execute(text(sql_query)), columns=columns) 
    return sources_df

def main():
  import argparse
  import os

  parser = argparse.ArgumentParser()
  
  parser.add_argument('-u', '--user', required=True, help='User for IEDB MySQL connection.')
  parser.add_argument('-p', '--password', required=True, help='User for IEDB MySQL connection.')
  parser.add_argument('-t', '--taxon_id', required=True, help='Taxon ID for the species to pull data for.')
  
  args = parser.parse_args()
  user = args.user
  password = args.password
  taxon_id = args.taxon_id

  # read in IEDB species data
  species_df = pd.read_csv('species.csv')
  species_id_to_name_map = dict(zip(species_df['Taxon ID'].astype(str), species_df['Species Label']))
  all_taxa_map = dict(zip(species_df['Taxon ID'].astype(str), species_df['All Taxa']))

  # get epitopes and source antigens
  Fetcher = DataFetcher(user, password)
  epitopes_df = Fetcher.get_epitopes(all_taxa_map[taxon_id])
  sources_df = Fetcher.get_sources(all_taxa_map[taxon_id])

  # create directory for species and taxon ID
  species_path = f'species/{taxon_id}-{species_id_to_name_map[taxon_id].replace(" ", "_")}'
  os.makedirs(species_path, exist_ok=True)

  # write epitopes and source antigens to files
  epitopes_df.to_csv(f'{species_path}/epitopes.csv', index=False)
  sources_df.to_csv(f'{species_path}/sources.csv', index=False)

if __name__ == '__main__':
  main()
