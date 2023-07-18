#!/usr/bin/env python3

import os
import pandas as pd
from pathlib import Path
from sqlalchemy import text

from sql_engine import create_sql_engine


class DataFetcher:
  def __init__(self, all_taxa: str) -> None:
    self.sql_engine = create_sql_engine()
    self.all_taxa = all_taxa.replace(';', ',')


  def get_data(self) -> tuple:
    """Get all epitopes and source antigens for a species together."""
    epitopes_df = self._get_epitopes()
    sources_df = self._get_sources()

    # keep only epitopes whose source antigen is in sources_df
    epitopes_df = epitopes_df[
      epitopes_df['Source Accession'].isin(sources_df['Accession'])
    ]

    return epitopes_df, sources_df


  def _get_epitopes(self) -> pd.DataFrame:
    """Get all epitopes for a species."""
    sql_query1 = f"""
                  SELECT object.mol1_seq, object.region, object.mol2_name, 
                         object.mol2_accession
                  FROM epitope, object
                  WHERE epitope.e_object_id = object.object_id
                  AND object.object_sub_type IN 
                      ("Peptide from protein", "Discontinuous protein residues")
                  AND object.organism2_id IN ({self.all_taxa});
                  """
    sql_query2 = f"""
                  SELECT object.mol1_seq, object.region, object.mol2_name, 
                         object.mol2_accession
                  FROM epitope, object
                  WHERE epitope.related_object_id = object.object_id
                  AND object.object_sub_type IN 
                      ("Peptide from protein", "Discontinuous protein residues")
                  AND object.organism2_id IN ({self.all_taxa});
                  """
    
    with self.sql_engine.connect() as connection:
      result1 = connection.execute(text(sql_query1))
      result2 = connection.execute(text(sql_query2))

      columns=[
        'Linear Sequence', 'Discontinuous Sequence', 'Source Name', 'Source Accession'
      ]
      df1 = pd.DataFrame(result1.fetchall(), columns=columns)
      df2 = pd.DataFrame(result2.fetchall(), columns=columns)
      
    epitopes_df = pd.concat([df1, df2], ignore_index=True)
    
    # combine linear and discontinuous sequences
    epitopes_df['Sequence'] = epitopes_df['Linear Sequence'].fillna(
      epitopes_df['Discontinuous Sequence']
    )
    epitopes_df.drop( # remove original sequence
      columns=['Linear Sequence', 'Discontinuous Sequence'], 
      inplace=True
    )
    # remove epitopes with no sequence
    epitopes_df.dropna(subset=['Sequence'], inplace=True)
    epitopes_df.drop_duplicates(inplace=True)
    
    return epitopes_df[['Sequence', 'Source Name', 'Source Accession']]

  def _get_sources(self) -> pd.DataFrame:
    """Get all source antigens for a species."""
    sql_query = f'SELECT accession, name, sequence '\
                f'FROM source WHERE organism_id IN ({self.all_taxa});'

    with self.sql_engine.connect() as connection:
      result = connection.execute(text(sql_query))
      sources_df = pd.DataFrame(result.fetchall(), columns=['Accession', 'Name', 'Sequence'])

    # remove sources with no sequence
    sources_df.dropna(subset=['Sequence'], inplace=True) 
    sources_df.drop_duplicates(inplace=True)
    sources_df['Length'] = sources_df['Sequence'].str.len()
    
    return sources_df
    

def run(taxon_id: int, all_taxa: str, species_name: str) -> None:
  """Run the data fetching pipeline for a species.

  Args:
    taxon_id: Taxon ID for species to pull data for.
    all_taxa: All children taxa for that species.
    species_name: Name of the species.
  """
  epitopes_df, sources_df = DataFetcher(all_taxa).get_data()

  # create directory for species and taxon ID
  species_path = Path(__file__).parent.parent / 'species' / f'{taxon_id}-{species_name.replace(" ", "_")}'
  os.makedirs(species_path, exist_ok=True)

  # write epitopes and source antigens to files
  epitopes_df.to_csv(species_path / 'epitopes.csv', index=False)
  sources_df.to_csv(species_path / 'sources.csv', index=False)


def main():
  import argparse
  
  parser = argparse.ArgumentParser()
  parser.add_argument(
    '-a', '--all_species',
    action='store_true',
    help='Build protein tree for all IEDB species.')
  parser.add_argument(
    '-t', '--taxon_id',
    type=int, 
    help='Taxon ID for species to pull data for.')
  
  args = parser.parse_args()

  species_df = pd.read_csv(Path(__file__).parent.parent / 'species.csv')
  valid_taxon_ids = species_df['Species Taxon ID'].tolist()
  
  # create maps for taxon ID to species name and all taxa
  species_name_map = dict(
    zip(
      species_df['Species Taxon ID'], 
      species_df['Species Name']
    )
  )
  all_taxa_map = dict(
    zip(
      species_df['Species Taxon ID'],
      species_df['All Taxa']
    )
  )  

  if not args.all_species and not args.taxon_id:
    print("Error: Please provide either --all_species (-a) or --taxon_id (-t).")
    return

  if args.all_species:
    for taxon_id in valid_taxon_ids:
      run(taxon_id, all_taxa_map[taxon_id], species_name_map[taxon_id])

  else:
    taxon_id = args.taxon_id
    run(taxon_id, all_taxa_map[taxon_id], species_name_map[taxon_id])

if __name__ == '__main__':
  main()
