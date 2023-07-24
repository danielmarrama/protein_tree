import pytest
import pandas as pd
from pathlib import Path

import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from select_proteome import ProteomeSelector
from assign_gene_protein import GeneAndProteinAssigner


data_path = Path(__file__).parent / "data"


@pytest.fixture(
  params=[
    [24, 'Shewanella putrefaciens', False, 'UP000008209'],
    [10042, 'Peromyscus maniculatus', True, 'UP000504601'],
    [334205, 'Nupapillomavirus 1', False, 'UP000006367']
  ]
)
def organism(request):
    return request.param


@pytest.fixture
def epitopes(organism) -> Path:
  taxon_id, species_name, _, _ = organism
  return data_path / f"{taxon_id}-{species_name.replace(' ', '_')}" / "epitopes.csv"


@pytest.fixture
def sources(organism) -> Path:
  taxon_id, species_name, _, _= organism
  return data_path / f"{taxon_id}-{species_name.replace(' ', '_')}" / "sources.csv"


def test_select_proteome(epitopes, organism):
  taxon_id, species_name, _, proteome_id = organism 

  epitopes_df = pd.read_csv(epitopes)
  selector = ProteomeSelector(taxon_id, species_name)
  proteome_data = selector.select_best_proteome(epitopes_df)
  assert proteome_data[0] == proteome_id


def test_assignments(epitopes, sources, organism):
  taxon_id, species_name, is_vertebrate, _ = organism
  
  epitopes_df = pd.read_csv(epitopes)
  sources_df = pd.read_csv(sources)

  Assigner = GeneAndProteinAssigner(taxon_id, species_name, is_vertebrate)
  assigner_data = Assigner.assign(sources_df, epitopes_df)
  
  assert assigner_data[2] == 1
  assert assigner_data[3] == 1