import unittest
import os
import pandas as pd

from pathlib import Path
from protein_tree.select_proteome import ProteomeSelector

class TestProteomeSelection(unittest.TestCase):
  def test_proteome_list_pull(self):
    selector = ProteomeSelector(taxon_id='12345', build_path=Path(__file__).parent / 'build')
    self.assertIsInstance(selector.proteome_list, pd.DataFrame)

  def test_gp_proteome_pull(self):
    selector = ProteomeSelector(taxon_id='1235689', build_path=Path(__file__).parent / 'build')
    print(selector.species_df.columns)
    selector._get_gp_proteome_to_fasta('UP000000216', '1235689')
    
    file_path = selector.build_path / 'species' / '12345-UP000000216' / '1235689_proteome.fasta'
    self.assertTrue(os.path.exists(file_path))

  def test_rest_proteome_pull(self):
    pass

  def test_taxonomy_protein_map_pull(self):
    pass

