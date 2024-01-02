import unittest
import pandas as pd

from pathlib import Path
from protein_tree.select_proteome import ProteomeSelector

class TestProteomeSelection(unittest.TestCase):
  def test_proteome_list_pull(self):
    selector = ProteomeSelector(taxon_id='12345', build_path=Path(__file__).parent / 'build')
    self.assertIsInstance(selector.proteome_list, pd.DataFrame)

  def test_gp_proteome_pull(self):
    pass

  def test_rest_proteome_pull(self):
    pass

  def test_taxonomy_protein_map_pull(self):
    pass

