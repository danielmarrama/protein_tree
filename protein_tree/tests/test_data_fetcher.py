from protein_tree.get_data import DataFetcher


def test_get_epitopes():
  fetcher = DataFetcher()
  result = fetcher.get_epitopes("319224") # Shewanella putrefaciens
  assert 'THEILWPSF' in result['Sequence'].values


def test_get_sources():
  fetcher = DataFetcher()
  result = fetcher.get_sources("319224") # Shewanella putrefaciens
  assert 'ABP77492.1' in result['Accession'].values