#!/usr/bin/env python3

import requests
import pandas as pd
from pathlib import Path
from sqlalchemy import text

class DataFetcher:
  def __init__(self, build_path: Path = Path(__file__).parent.parent / 'build') -> None:
    self.build_path = build_path
    
  def get_all_data(self) -> None:
    """Get both peptides and peptide sources tables from the IEDB backend. Also, get
    allergen data from the IUIS allergen nomenclature database.
    """

    from sql_engine import create_sql_engine
    self.sql_engine = create_sql_engine()
    peptides_df = self._get_peptide_table()
    sources_df = self._get_source_table()

    # keep only peptides whose source is in sources_df
    peptides_df = peptides_df[
      peptides_df['Source Accession'].isin(sources_df['Accession'])
    ]

    # get allergen data
    url = 'http://www.allergen.org/csv.php?table=joint'
    allergen_df = pd.read_csv(url)
    
    # write data to TSV files
    peptides_df.to_csv(self.build_path / 'iedb' / 'peptides.tsv', sep='\t', index=False)
    sources_df.to_csv(self.build_path / 'iedb' / 'sources.tsv', sep='\t', index=False)
    allergen_df.to_csv(self.build_path / 'arborist' / 'allergens.tsv', sep='\t', index=False)

  def _get_peptide_table(self) -> pd.DataFrame:
    """Get peptides from peptide/object tables from the IEDB backend."""
    sql_query1 = f"""
      SELECT object.mol1_seq, object.region, object.mol2_name, 
             object.mol2_accession, object.organism2_id
      FROM peptide, object
      WHERE peptide.e_object_id = object.object_id
      AND object.object_sub_type IN 
        ("Peptide from protein", "Discontinuous protein residues")
      """
    sql_query2 = f"""
      SELECT object.mol1_seq, object.region, object.mol2_name,
             object.mol2_accession, object.organism2_id
      FROM peptide, object
      WHERE peptide.related_object_id = object.object_id
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
      
    peptides_df = pd.concat([df1, df2], ignore_index=True)
    
    # combine linear and discontinuous sequences
    peptides_df['Sequence'] = peptides_df['Linear Sequence'].fillna(
      peptides_df['Discontinuous Sequence']
    )
    peptides_df.drop( # remove original sequence
      columns=['Linear Sequence', 'Discontinuous Sequence'], 
      inplace=True
    )
    # remove peptides with no sequence
    peptides_df.dropna(subset=['Sequence'], inplace=True)

    # limit to peptides with 5 or more amino acids
    peptides_df = peptides_df[peptides_df['Sequence'].str.len() >= 5]

    # make peptide sequences uppercase
    peptides_df['Sequence'] = peptides_df['Sequence'].str.upper()

    peptides_df.drop_duplicates(inplace=True)
  
    return peptides_df[['Sequence', 'Source Name', 'Source Accession', 'Organism ID']]

  def _get_source_table(self) -> pd.DataFrame:
    """Get all peptide sources from source table in IEDB backend."""
    sql_query = f'SELECT accession, name, sequence, organism_id FROM source;'

    with self.sql_engine.connect() as connection:
      result = connection.execute(text(sql_query))
      sources_df = pd.DataFrame(
        result.fetchall(), columns=['Accession', 'Name', 'Sequence', 'Organism ID']
      )

    # remove sources with no sequence
    # sources_df.dropna(subset=['Sequence'], inplace=True) 
    sources_df.drop_duplicates(inplace=True)
    sources_df['Length'] = sources_df['Sequence'].str.len()
    
    return sources_df
  
  def get_all_peptides(self) -> pd.DataFrame:
    """Get all peptides from the written file."""
    return pd.read_csv(self.build_path / 'iedb' / 'peptide.tsv', sep='\t')
  
  def get_all_sources(self) -> pd.DataFrame:
    """Get all peptide sources from the written file."""
    return pd.read_csv(self.build_path / 'iedb' / 'peptide_source.tsv', sep='\t')

  def get_peptides_for_species(
    self, all_peptides: pd.DataFrame, all_taxa: list
  ) -> pd.DataFrame:
    """Get peptides from the written file only for a specific species.
    
    Args:
      all_peptides: list of all peptides from the backend. 
      all_taxa: list of all active children taxa for a species.
    """
    return all_peptides[all_peptides['Organism ID'].isin(all_taxa)]

  def get_sources_for_species(
    self, all_sources: pd.DataFrame, accessions: list
  ) -> pd.DataFrame:
    """Get peptide sources from the written file only for a specific species.
    
    Args:
      all_sources: list of all peptide sources from the backend.
      all_taxa: list of all active children taxa for a species.
    """
    return all_sources[all_sources['Accession'].isin(accessions)]

  @staticmethod
  def update_species():
    """Update species table "active-species-tsv" fron OntoDev site."""
    url = 'https://nb.ontodev.com/active_species.tsv?limit=3000&offset=0'
    
    r = requests.get(url)
    r.raise_for_status()

    file_path = Path(__file__).parent.parent / 'data' / 'active-species.tsv'
    with open(file_path, 'wb') as f:
      f.write(r.content)

if __name__ == '__main__':
  import argparse
  
  parser = argparse.ArgumentParser()
  parser.add_argument(
    '-t', '--taxon_id',
    type=int, 
    help='Write separate files for a species with its taxon ID.'
  )
  parser.add_argument(
    '-d', '--data_path',
    type=str,
    default=Path(__file__).parent.parent / 'data',
    help='Directory to write data to.'
  )

  args = parser.parse_args()
  data_path = Path(args.data_path)
  
  species_df = pd.read_csv(data_path / 'active-species.tsv', sep='\t')
  all_taxon_ids = species_df['Species ID'].tolist()
  
  # create maps for taxon ID to species name and all taxa
  species_name_map = dict(
    zip(
      species_df['Species ID'], 
      species_df['Species Label']
    )
  )
  all_taxa_map = dict(
    zip(
      species_df['Species ID'],
      species_df['Active Taxa']
    )
  )  

  if args.taxon_id:
    files_exist = (
      (data_path / 'iedb' / 'peptides.tsv').exists() and
      (data_path / 'iedb' / 'sources.tsv').exists() and
      (data_path / 'arborist' / 'allergens.tsv').exists()
    )
    if not files_exist:
      print('Getting all data...')
      DataFetcher().get_all_data()
      print('All data written.')

    all_peptides = DataFetcher().get_all_peptides()
    all_sources = DataFetcher().get_all_sources()

    taxon_id = args.taxon_id
    assert taxon_id in all_taxon_ids, f'{taxon_id} is not a valid taxon ID.'
    
    species_name = species_name_map[taxon_id]

    species_path = data_path / 'species' / f'{taxon_id}' # directory to write species data
    species_path.mkdir(parents=True, exist_ok=True)

    print(f'Writing separate files for {species_name}...')
    all_taxa = [int(taxon) for taxon in all_taxa_map[taxon_id].split(', ')]
    
    peptides_df = DataFetcher().get_peptides_for_species(all_peptides, all_taxa)
    sources_df = DataFetcher().get_sources_for_species(
      all_sources, peptides_df['Source Accession'].tolist()
    )
    
    peptides_df.to_csv(species_path / 'peptides.tsv', sep='\t', index=False)
    sources_df.to_csv(species_path / 'sources.tsv', sep='\t', index=False)
    print(f'peptides and sources for {species_name} written.')

  else:
    print('Getting all data...')
    DataFetcher().get_all_data()
    print('All data written.')
