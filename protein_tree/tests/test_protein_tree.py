import pytest
import pandas as pd
from pathlib import Path

import os
import sys
import glob

# add path to parent directory to sys.path so that we can import modules
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from select_proteome import ProteomeSelector
from assign_gene_protein import GeneAndProteinAssigner


data_path = Path(__file__).parent / "data"


@pytest.fixture(
  params=[
    [106, 'Runella slithyformis', False, 'UP000000493'],
    [10042, 'Peromyscus maniculatus', True, 'UP000504601'],
    [334205, 'Nupapillomavirus 1', False, 'UP000006367']
  ]
)
def organism(request):
    return request.param


@pytest.fixture
def epitopes(organism) -> Path:
  taxon_id, species_name, _, _ = organism
  return data_path / 'species' / f"{taxon_id}-{species_name.replace(' ', '_')}" / "epitopes.csv"


@pytest.fixture
def sources(organism) -> Path:
  taxon_id, species_name, _, _= organism
  return data_path / 'species' / f"{taxon_id}-{species_name.replace(' ', '_')}" / "sources.csv"


@pytest.fixture(scope='session', autouse=True)
def cleanup_after_all_tests():
  yield
  for file in glob.glob(str(data_path / 'species' / '**' / '*proteome*'), recursive=True):
    os.remove(file)


def test_select_proteome(epitopes, organism):
  taxon_id, species_name, _, proteome_id = organism 

  epitopes_df = pd.read_csv(epitopes)
  Selector = ProteomeSelector(taxon_id, species_name, data_path)
  proteome_data = Selector.select_best_proteome(epitopes_df)
  Selector.proteome_to_csv()

  assert proteome_data[0] == proteome_id


def test_assignments(epitopes, sources, organism):
  taxon_id, species_name, is_vertebrate, _ = organism
  
  epitopes_df = pd.read_csv(epitopes)
  sources_df = pd.read_csv(sources)

  Assigner = GeneAndProteinAssigner(taxon_id, species_name, is_vertebrate, data_path)
  _, epitope_assignments, source_assignments = Assigner.assign(sources_df, epitopes_df)
  
  epitopes_expected = pd.read_csv(data_path / 'species' / f"{taxon_id}-{species_name.replace(' ', '_')}" / 'epitope_assignments.csv')
  sources_expected = pd.read_csv(data_path / 'species' / f"{taxon_id}-{species_name.replace(' ', '_')}" / 'source_assignments.csv')

  assert epitopes_expected.equals(epitope_assignments)
  assert sources_expected.equals(source_assignments)
