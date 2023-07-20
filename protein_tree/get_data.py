#!/usr/bin/env python3

import os
import pandas as pd
from pathlib import Path
from sqlalchemy import text

from sql_engine import create_sql_engine


class DataFetcher:
  def __init__(self) -> None:
    self.data_dir = Path(__file__).parent.parent / 'data'
    self.sql_engine = create_sql_engine()

    self.epitopes_df = pd.read_csv(self.data_dir / 'epitopes.csv')
    self.sources_df = pd.read_csv(self.data_dir / 'sources.csv')


  def get_all_data(self) -> None:
    """Get all epitopes and source antigens tables."""
    epitopes_df = self._get_epitope_table()
    sources_df = self._get_source_table()

    # keep only epitopes whose source antigen is in sources_df
    epitopes_df = epitopes_df[
      epitopes_df['Source Accession'].isin(sources_df['Accession'])
    ]

    # get allergen data
    url = 'http://www.allergen.org/csv.php?table=joint'
    allergen_df = pd.read_csv(url)
    
    # write tables to fiels
    epitopes_df.to_csv(self.data_dir / 'epitopes.csv', index=False)
    sources_df.to_csv(self.data_dir / 'sources.csv', index=False)
    allergen_df.to_csv(self.data_dir / 'allergens.csv', index=False)


  def _get_epitope_table(self) -> pd.DataFrame:
    """Get epitopes table from the IEDB backend."""
    sql_query1 = f"""
      SELECT object.mol1_seq, object.region, object.mol2_name, 
             object.mol2_accession, object.organism2_id
      FROM epitope, object
      WHERE epitope.e_object_id = object.object_id
      AND object.object_sub_type IN 
          ("Peptide from protein", "Discontinuous protein residues")
      """
    sql_query2 = f"""
      SELECT object.mol1_seq, object.region, object.mol2_name,
             object.mol2_accession, object.organism2_id
      FROM epitope, object
      WHERE epitope.related_object_id = object.object_id
      AND object.object_sub_type IN 
          ("Peptide from protein", "Discontinuous protein residues")
      """
    
    with self.sql_engine.connect() as connection:
      result1 = connection.execute(text(sql_query1))
      result2 = connection.execute(text(sql_query2))

      columns=[
        'Linear Sequence', 'Discontinuous Sequence', 
        'Source Name', 'Source Accession', 'Organism ID'
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
    
    return epitopes_df[['Sequence', 'Source Name', 'Source Accession', 'Organism ID']]


  def _get_source_table(self) -> pd.DataFrame:
    """Get all source antigens for a species."""
    sql_query = f'SELECT accession, name, sequence, organism_id FROM source;'

    with self.sql_engine.connect() as connection:
      result = connection.execute(text(sql_query))
      sources_df = pd.DataFrame(
        result.fetchall(), columns=['Accession', 'Name', 'Sequence', 'Organism ID']
      )

    # remove sources with no sequence
    sources_df.dropna(subset=['Sequence'], inplace=True) 
    sources_df.drop_duplicates(inplace=True)
    sources_df['Length'] = sources_df['Sequence'].str.len()
    
    return sources_df
  

  def get_epitopes_for_species(self, all_taxa: list) -> pd.DataFrame:
    """Get epitopes for a species.
    
    Args:
      all_taxa: list of all active children taxa for a species.
    """
    return self.epitopes_df[self.epitopes_df['Organism ID'].isin(all_taxa)]


  def get_sources_for_species(self, all_taxa: list) -> pd.DataFrame:
    """Get source antigens for a species.
    
    Args:
      all_taxa: list of all active children taxa for a species.
    """
    return self.sources_df[self.sources_df['Organism ID'].isin(all_taxa)]


def main():
  import argparse
  
  parser = argparse.ArgumentParser()
  parser.add_argument(
    '-t', '--taxon_id',
    type=int, 
    help='Write separate files for a species with its taxon ID.'
  )
  
  args = parser.parse_args()

  data_dir = Path(__file__).parent.parent / 'data'
  species_df = pd.read_csv(data_dir / 'species.csv')
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

  if args.taxon_id:
    if not (data_dir / 'epitopes.csv').exists() or not (data_dir / 'sources.csv').exists():
      print('Getting all data...')
      DataFetcher().get_all_data()
      print('All data written.')

    taxon_id = args.taxon_id
    assert taxon_id in valid_taxon_ids, f'{taxon_id} is not a valid taxon ID.'

    species_name = species_name_map[taxon_id]
    species_dir = data_dir / 'species' / f'{taxon_id}-{species_name.replace(" ", "_")}'
    species_dir.mkdir(parents=True, exist_ok=True)

    print(f'Writing separate files for {species_name}...')
    all_taxa = [int(taxon) for taxon in all_taxa_map[taxon_id].split(';')]
    epitopes_df = DataFetcher().get_epitopes_for_species(all_taxa)
    sources_df = DataFetcher().get_sources_for_species(all_taxa)
    
    epitopes_df.to_csv(species_dir / 'epitopes.csv', index=False)
    sources_df.to_csv(species_dir / 'sources.csv', index=False)
    print(f'Epitopes and sources for {species_name} written.')

  else:
    print('Getting all data...')
    DataFetcher().get_all_data()
    print('All data written.')


if __name__ == '__main__':
  main()
