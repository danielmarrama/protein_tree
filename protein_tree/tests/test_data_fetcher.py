from unittest.mock import patch
from protein_tree.get_data import DataFetcher
import pandas as pd

# Create some dummy data to represent the results from the database query
dummy_epitopes = pd.DataFrame({
    'Sequence': ['THEILWPSF', 'OTHERSEQUENCE'],
    'Source Name': ['Source 1', 'Source 2'],
    'Source Accession': ['ACC1', 'ACC2'],
})

dummy_sources = pd.DataFrame({
    'Accession': ['ABP77492.1', 'OTHERACC'],
    'Name': ['Name 1', 'Name 2'],
    'Sequence': ['SEQUENCE1', 'SEQUENCE2'],
})

# Patch the get_epitopes and get_sources methods in the DataFetcher class
@patch.object(DataFetcher, 'get_epitopes', return_value=dummy_epitopes)
@patch.object(DataFetcher, 'get_sources', return_value=dummy_sources)
def test_get_data(mock_get_epitopes, mock_get_sources):
  fetcher = DataFetcher()
  result = fetcher.get_epitopes("319224")  # Shewanella putrefaciens
  assert 'THEILWPSF' in result['Sequence'].values

  result = fetcher.get_sources("319224")  # Shewanella putrefaciens
  assert 'ABP77492.1' in result['Accession'].values
