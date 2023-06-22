import pytest
import pandas as pd
from pathlib import Path

from protein_tree.select_proteome import ProteomeSelector


@pytest.fixture
def epitopes() -> Path:
  return Path(__file__).parent / "data" / "epitopes.csv"


def test_select_proteome(epitopes):
  epitopes_df = pd.read_csv(epitopes)
  selector = ProteomeSelector('24', 'Shewanella putrefaciens')
  proteome_data = selector.select_best_proteome(epitopes_df)
  assert proteome_data[0] == 'UP000008209'
