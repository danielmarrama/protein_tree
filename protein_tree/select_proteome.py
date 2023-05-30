#!/usr/bin/env python3

import re
import os
import requests
import pandas as pd
from pathlib import Path

from pepmatch import Preprocessor, Matcher


class ProteomeSelector:
  def __init__(self, taxon_id, species_name):
    self.taxon_id = taxon_id

    parent_dir = Path(__file__).resolve().parent.parent
    species_dir = f'{taxon_id}-{species_name.replace(" ", "_")}'
    self.species_path = parent_dir / 'species' / species_dir

    self.species_df = pd.read_csv(f'{parent_dir}/species.csv') # IEBD species data
    self.metrics_df = pd.read_csv(f'{parent_dir}/metrics.csv') # protein tree metrics

    self.proteome_list = self._get_proteome_list() # get all candidate proteomes
    self.num_of_proteomes = len(self.proteome_list) + 1 # +1 for all proteins option

  def select_best_proteome(self, epitopes_df: pd.DataFrame) -> list:
    """
    Select the best proteome to use for a species. Return the proteome ID, 
    proteome taxon, and proteome type.
    
    Check UniProt for all candidate proteomes associated with that
    taxon. Then do the following checks:

    1. Are there any representative proteomes?
    2. Are there any reference proteomes?
    3. Are there any non-redudant proteomes?
    4. Are there any other proteomes?

    If yes to any of the above, check if there are multiples and do epitope
    search for tie breaks.

    If no to all of the above, then get every protein associated with
    the taxon ID using the get_all_proteins method.

    Args:
      epitopes_df (pd.DataFrame): DataFrame of epitopes for the species to use for tie breaks.
    """
    if (self.species_path / 'proteome.fasta').exists():
      idx = self.metrics_df['Species Taxon ID'] == self.taxon_id
      proteome_id = self.metrics_df[idx]['Proteome ID'].iloc[0]
      proteome_taxon = self.metrics_df[idx]['Proteome Taxon'].iloc[0]
      proteome_type = self.metrics_df[idx]['Proteome Type'].iloc[0]
      
      return proteome_id, proteome_taxon, proteome_type
    
    else:
      self.species_path.mkdir(parents=True, exist_ok=True)

    if self.proteome_list.empty:
      print('No proteomes found. Fetching orphan proteins.')
      self._get_all_proteins()
      return 'None', self.taxon_id, 'All-proteins'

    if self.proteome_list['isRepresentativeProteome'].any():
      print('Found representative proteome(s).')
      proteome_type = 'Representative'
      self.proteome_list = self.proteome_list[self.proteome_list['isRepresentativeProteome']]
      proteome_id, proteome_taxon = self._get_proteome_with_most_matches(epitopes_df)
      self._get_gp_proteome_to_fasta(proteome_id, proteome_taxon)
    
    elif self.proteome_list['isReferenceProteome'].any():
      print('Found reference proteome(s).')
      proteome_type = 'Reference'
      self.proteome_list = self.proteome_list[self.proteome_list['isReferenceProteome']]
      proteome_id, proteome_taxon = self._get_proteome_with_most_matches(epitopes_df)
      self._get_gp_proteome_to_fasta(proteome_id, proteome_taxon)

    elif 'redundantTo' not in self.proteome_list.columns:
      print('Found other proteome(s).')
      proteome_type = 'Other'
      proteome_id, proteome_taxon = self._get_proteome_with_most_matches(epitopes_df)
    
    elif self.proteome_list['redundantTo'].isna().any():
      print('Found non-redundant proteome(s).')
      proteome_type = 'Non-redundant'
      self.proteome_list = self.proteome_list[self.proteome_list['redundantTo'].isna()]
      proteome_id, proteome_taxon = self._get_proteome_with_most_matches(epitopes_df)
    
    else:
      print('Found other proteome(s).')
      proteome_type = 'Other'
      proteome_id, proteome_taxon = self._get_proteome_with_most_matches(epitopes_df)

    self._remove_other_proteomes(proteome_id)
    
    # sanity check to make sure proteome.fasta is not empty
    if (self.species_path / 'proteome.fasta').stat().st_size == 0:
      proteome_id = 'None'
      proteome_taxon = self.taxon_id
      proteome_type = 'All-proteins'
      self._get_all_proteins() # get all orphan proteins if proteome.fasta is empty

    return proteome_id, proteome_taxon, proteome_type

  def proteome_to_csv(self) -> None:
    """
    Write the proteome data for a species to a CSV file for later use.
    """
    from Bio import SeqIO
    
    if not (self.species_path / 'proteome.fasta').exists():
      return

    regexes = {
      'protein_id': re.compile(r"\|([^|]*)\|"),     # between | and |
      'protein_name': re.compile(r"\s(.+?)\sOS"),   # between space and space before OS
      'gene': re.compile(r"GN=(.+?)\s"),            # between GN= and space
      'pe_level': re.compile(r"PE=(.+?)\s"),        # between PE= and space
    }

    proteins = list(SeqIO.parse(f'{self.species_path}/proteome.fasta', 'fasta'))
    gp_proteome_path = self.species_path / 'gp_proteome.fasta'
    if gp_proteome_path.exists():
      gp_ids = [str(protein.id.split('|')[1]) for protein in list(SeqIO.parse(gp_proteome_path, 'fasta'))]
    else:
      gp_ids = []

    proteome_data = []
    for protein in proteins:
      metadata = []
      for key in regexes:
        match = regexes[key].search(str(protein.description))

        if match:
          metadata.append(match.group(1))
        else:
          if key == 'protein_id':
            metadata.append(str(protein.id))
          elif key == 'pe_level':
            metadata.append(0)
          else:
            metadata.append('')
      
      gp = 1 if protein.id.split('|')[1] in gp_ids else 0
      metadata.append(gp)
      metadata.append(str(protein.seq))
      metadata.append(protein.id.split('|')[0])
      
      proteome_data.append(metadata)
    
    columns = ['Protein ID', 'Protein Name', 'Gene', 'Protein Existence Level', 'Gene Priority', 'Sequence', 'Database']
    proteome = pd.DataFrame(proteome_data, columns=columns)
    proteome = proteome[['Database', 'Gene', 'Protein ID', 'Protein Name', 'Protein Existence Level', 'Gene Priority', 'Sequence']]
    proteome.to_csv(f'{self.species_path}/proteome.csv', index=False)


  def _get_proteome_list(self) -> pd.DataFrame:
    """
    Get a list of proteomes for a species from the UniProt API.
    Check for proteome_type:1 first, which are the representative or
    reference proteomes.

    If there are no proteomes, return empty DataFrame.
    """
    # URL to get proteome list for a species - use proteome_type:1 first
    url = f'https://rest.uniprot.org/proteomes/stream?format=xml&query=(proteome_type:1)AND(taxonomy_id:{self.taxon_id})'
    
    try:
      proteome_list = pd.read_xml(requests.get(url).text)
    except ValueError:
      try: # delete proteome_type:1 from URL and try again
        url = url.replace('(proteome_type:1)AND', '')
        proteome_list = pd.read_xml(requests.get(url).text)
      except ValueError: # if there are no proteomes, return empty DataFrame
        return pd.DataFrame()

    # remove the namespace from the columns
    proteome_list.columns = [x.replace('{http://uniprot.org/proteome}', '') for x in proteome_list.columns]
    return proteome_list

  def _get_all_proteins(self) -> None:
    """
    Get every protein associated with a taxon ID on UniProt.
    Species on UniProt will have a proteome, but not every protein is
    stored within those proteomes. There is a way to get every protein
    using the taxonomy part of UniProt. 
    """
    # URL link to all proteins for a species - size = 500 proteins at a time
    url = f'https://rest.uniprot.org/uniprotkb/search?format=fasta&query=taxonomy_id:{self.taxon_id}&size=500'

    # loop through all protein batches and write proteins to FASTA file
    for batch in self._get_protein_batches(url):
      with open(f'{self.species_path}/proteome.fasta', 'a') as f:
        f.write(batch.text)

  def _get_protein_batches(self, batch_url: str) -> requests.Response:
    """
    Get a batch of proteins from UniProt API because it limits the
    number of proteins you can get at once. Yield each batch until the 
    URL link is empty.
    
    Args:
      batch_url (str): URL to get all proteins for a species.
    """
    while batch_url:
      r = requests.get(batch_url)
      r.raise_for_status()
      yield r
      batch_url = self._get_next_link(r.headers)

  def _get_next_link(self, headers: dict) -> str:
    """
    UniProt will provide a link to the next batch of proteins when getting
    all proteins for a species' taxon ID.
    We can use a regular expression to extract the URL from the header.

    Args:
      headers (dict): Headers from UniProt API response.
    """
    re_next_link = re.compile(r'<(.+)>; rel="next"') # regex to extract URL
    if 'Link' in headers:
      match = re_next_link.match(headers['Link'])
      if match:
        return match.group(1)

  def _get_gp_proteome_to_fasta(self, proteome_id: str, proteome_taxon: str) -> None:
    """
    Write the gene priority proteome to a file. 
    This is only for representative and reference proteomes.
    Depending on the species group, the FTP URL will be different.

    Args:
      proteome_id (str): Proteome ID.
      proteome_taxon (str): Taxon ID for the proteome.
    """
    import gzip
    
    idx = self.species_df['Species Taxon ID'] == self.taxon_id
    group = self.species_df[idx]['Group'].iloc[0]
    ftp_url = f'https://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/reference_proteomes/'
    
    if group == 'Archaea':
      ftp_url += f'Archaea/{proteome_id}/{proteome_id}_{proteome_taxon}.fasta.gz'
    elif group == 'Bacteria':
      ftp_url += f'Bacteria/{proteome_id}/{proteome_id}_{proteome_taxon}.fasta.gz'
    elif group == 'Eukaryota':
      ftp_url += f'Eukaryota/{proteome_id}/{proteome_id}_{proteome_taxon}.fasta.gz'
    elif group == 'Viruses':
      ftp_url += f'Viruses/{proteome_id}/{proteome_id}_{proteome_taxon}.fasta.gz'
    
    r = requests.get(ftp_url, stream=True)
    try:
      r.raise_for_status()
    except:
      return

    # unzip the request and write the gene priority proteome to a file
    with open(f'{self.species_path}/gp_proteome.fasta', 'wb') as f:
      f.write(gzip.open(r.raw, 'rb').read())

  def _get_proteome_to_fasta(self, proteome_id: str) -> None:
    """
    Get the FASTA file for a proteome from UniProt API.
    Include all isoforms and do not compress the file.

    Args:
      proteome_id (str): UniProt Proteome ID.
    """
    url = f'https://rest.uniprot.org/uniprotkb/stream?compressed=false&format=fasta&includeIsoform=true&query=(proteome:{proteome_id})'
    r = requests.get(url)
    r.raise_for_status()
    with open(f'{self.species_path}/{proteome_id}.fasta', 'w') as f:
      f.write(r.text)

  def _get_proteome_with_most_matches(self, epitopes_df: pd.DataFrame) -> tuple:
    """
    Get the proteome ID and true taxon ID associated with
    that proteome with the most epitope matches in case there is a tie.
    We get the true taxon so we can extract the data from the FTP server
    if needed.

    Args:
      epitopes_df (pd.DataFrame): DataFrame of epitopes for the species.
    """
    if self.num_of_proteomes <= 2:
      self._get_proteome_to_fasta(self.proteome_list['upid'].iloc[0])
      return self.proteome_list['upid'].iloc[0], self.proteome_list['taxonomy'].iloc[0]

    epitopes_df = epitopes_df[epitopes_df['Sequence'].notna()] 
    
    # TODO: be able to check for discontinuous epitopes
    epitopes = [epitope for epitope in epitopes_df['Sequence'] if not any(char.isdigit() for char in epitope)]
    
    match_counts = {} # keep track of # of epitope matches for each proteome
    for proteome_id in list(self.proteome_list['upid']):
      self._get_proteome_to_fasta(proteome_id)
      
      Preprocessor(
        proteome = f'{self.species_path}/{proteome_id}.fasta',
        preprocessed_files_path = f'{self.species_path}',
      ).sql_proteome(k = 5)
      
      matches_df = Matcher(
        query = epitopes, 
        proteome_file = f'{self.species_path}/{proteome_id}.fasta', 
        max_mismatches = 0, 
        k = 5,
        preprocessed_files_path = f'{self.species_path}',
        output_format='dataframe'
      ).match()
      
      matches_df.drop_duplicates(subset=['Query Sequence'], inplace=True)
      
      try:
        match_counts[proteome_id] = matches_df['Matched Sequence'].dropna().count()
      except KeyError: # in case there are no matches
        match_counts[proteome_id] = 0 

    # select the proteome ID and proteome taxon with the most matches
    proteome_id = max(match_counts, key=match_counts.get)
    proteome_taxon = self.proteome_list[self.proteome_list['upid'] == proteome_id]['taxonomy'].iloc[0]

    return proteome_id, proteome_taxon

  def _remove_other_proteomes(self, proteome_id: str) -> None:
    """
    Remove the proteome FASTA files that are not the chosen proteome for that
    species. Also, remove the .db files and rename the chosen proteome to 
    "proteome.fasta".

    Args:
      proteome_id (str): Proteome ID of the chosen proteome.
    """
    proteome_list_to_remove = self.proteome_list[self.proteome_list['upid'] != proteome_id]
    for i in list(proteome_list_to_remove['upid']):
      os.remove(self.species_path / f'{i}.fasta')
      os.remove(self.species_path / f'{i}.db')
    
    # rename the chosen proteome to proteome.fasta and remove the .db file
    os.rename(f'{self.species_path}/{proteome_id}.fasta', f'{self.species_path}/proteome.fasta')
    if self.num_of_proteomes > 2: # there is only a .db file if there is more than one proteome
      os.remove(f'{self.species_path}/{proteome_id}.db')


def main():
  import argparse
  from get_data import DataFetcher

  parser = argparse.ArgumentParser()
  
  parser.add_argument(
    '-u', '--user', 
    required=True, 
    help='User for IEDB MySQL connection.')
  parser.add_argument(
    '-p', '--password', 
    required=True, 
    help='Password for IEDB MySQL connection.')
  parser.add_argument(
    '-a', '--all_species', 
    action='store_true', 
    help='Build protein tree for all IEDB species.')
  parser.add_argument(
    '-t', '--taxon_id', 
    required=True, 
    help='Taxon ID for the species to pull data for.')
  
  args = parser.parse_args()
  user = args.user
  password = args.password
  all_species = args.all_species
  taxon_id = args.taxon_id

  species_df = pd.read_csv('species.csv') # IEDB species data
  metrics_df = pd.read_csv('metrics.csv') # protein tree metrics data
  valid_taxon_ids = species_df['Species Taxon ID'].tolist()

  # taxa and name mapppings
  all_taxa_map = dict(zip(species_df['Species Taxon ID'].astype(str), species_df['All Taxa']))
  species_id_to_name_map = dict(zip(species_df['Species Taxon ID'].astype(str), species_df['Species Name']))

  if all_species:
    proteomes = {}
    for taxon_id in valid_taxon_ids:
      Fetcher = DataFetcher(user, password)
      epitopes_df = Fetcher.get_epitopes(all_taxa_map[taxon_id])

      print(f'Selecting best proteome for {species_id_to_name_map[taxon_id]} (Taxon ID: {taxon_id}).')
      Selector = ProteomeSelector(taxon_id, species_id_to_name_map[taxon_id], species_df, metrics_df)
      print(f'Number of candidate proteomes: {Selector.num_of_proteomes}')

      proteome_data = Selector.select_best_proteome(epitopes_df)
      Selector.proteome_to_csv()
      
      print(f'Proteome ID: {proteome_data[0]}')
      print(f'Proteome taxon: {proteome_data[1]}')
      print(f'Proteome type: {proteome_data[2]}')
    
      # update metrics data
      metrics_df.loc[metrics_df['Species Taxon ID'] == int(t_id), 'Proteome ID'] = proteome_data[0]
      metrics_df.loc[metrics_df['Species Taxon ID'] == int(t_id), 'Proteome Taxon'] = proteome_data[1]
      metrics_df.loc[metrics_df['Species Taxon ID'] == int(t_id), 'Proteome Type'] = proteome_data[2]

      metrics_df.to_csv('metrics.csv', index=False)

  else: # one species at a time
    assert taxon_id in valid_taxon_ids, f'{taxon_id} is not a valid taxon ID.'
    
    Fetcher = DataFetcher(user, password)
    epitopes_df = Fetcher.get_epitopes(all_taxa_map[taxon_id])

    print(f'Selecting best proteome for {species_id_to_name_map[taxon_id]} (Taxon ID: {taxon_id}).')
    Selector = ProteomeSelector(taxon_id, species_id_to_name_map[taxon_id], species_df, metrics_df)
    print(f'Number of candidate proteomes: {Selector.num_of_proteomes}')

    proteome_data = Selector.select_best_proteome(epitopes_df)
    Selector.proteome_to_csv()

    print(f'Proteome ID: {proteome_data[0]}')
    print(f'Proteome taxon: {int(proteome_data[1])}')
    print(f'Proteome type: {proteome_data[2]}')

    # update metrics data
    metrics_df.loc[metrics_df['Species Taxon ID'] == int(taxon_id), 'Proteome ID'] = proteome_data[0]
    metrics_df.loc[metrics_df['Species Taxon ID'] == int(taxon_id), 'Proteome Taxon'] = int(proteome_data[1])
    metrics_df.loc[metrics_df['Species Taxon ID'] == int(taxon_id), 'Proteome Type'] = proteome_data[2]

    metrics_df.to_csv('metrics.csv', index=False)

if __name__ == '__main__':
  main()