import pandas as pd
from unittest.mock import patch
from protein_tree.get_data import DataFetcher
from protein_tree.select_proteome import ProteomeSelector


# Create some dummy data to represent the results from the database query
dummy_epitopes = pd.DataFrame({
    'Sequence': ['THEILWPSF'],
    'Source Name': ['cytochrome c oxidase, subunit III'],
    'Source Accession': ['ABP77492.1'],
})

dummy_sources = pd.DataFrame({
    'Accession': ['ABP77492.1'],
    'Name': ['cytochrome c oxidase, subunit III'],
})

# Patch the get_epitopes and get_sources methods in the DataFetcher class
# VPN is required to access the database so we don't want to actually run
# the queries in the tests
@patch('protein_tree.get_data.DataFetcher')
def test_get_data(MockDataFetcher):
  mock_fetcher = MockDataFetcher.return_value
  mock_fetcher.get_epitopes.return_value = dummy_epitopes
  mock_fetcher.get_sources.return_value = dummy_sources

  epitopes = mock_fetcher.get_epitopes('319224')
  assert 'THEILWPSF' in epitopes['Sequence'].values

  sources = mock_fetcher.get_sources('319224')
  assert 'ABP77492.1' in sources['Accession'].values

@patch('protein_tree.get_data.DataFetcher')
def test_select_proteome(MockDataFetcher):
  mock_fetcher = MockDataFetcher.return_value
  mock_fetcher.get_epitopes.return_value = dummy_epitopes
  epitopes = mock_fetcher.get_epitopes('319224')

  selector = ProteomeSelector('24', 'Shewanella putrefaciens')
  proteome_data = selector.select_best_proteome(epitopes)
  assert proteome_data[0] == 'UP000008209'

