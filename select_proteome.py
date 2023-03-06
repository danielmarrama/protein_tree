#!/usr/bin/env python3

import warnings
warnings.filterwarnings('ignore')

import re
import os
import gzip
import argparse
import pandas as pd
import requests
from requests.adapters import HTTPAdapter, Retry
from pepmatch import Preprocessor, Matcher
from Bio import SeqIO

from get_data import DataFetcher


class ProteomeSelector:
  def __init__(self, taxon_id):
    self.species_df = pd.read_csv('species.csv')
    self.taxon_id = taxon_id
    self.proteome_list = self.get_proteome_list()
    self.num_of_proteomes = len(self.proteome_list)

  def get_proteome_list(self):
    """
    Get a list of proteomes for a species from the UniProt API.

    If there are no proteomes, return None.
    """
    # URL to get proteome list for a species - use proteome_type:1 first to get reference proteomes
    url = f'https://rest.uniprot.org/proteomes/stream?format=xml&query=(proteome_type:1)AND(taxonomy_id:{self.taxon_id})'
    
    # read in proteome list for the species from the UniProt API - try to get
    # the proteome list with proteome_type:1, which are the reference proteomes
    # if there are no reference proteomes, then try to get the other proteomes
    try:
      proteome_list = pd.read_xml(requests.get(url).text)
    except ValueError:
      try:
        url = url.replace('(proteome_type:1)AND', '')
        proteome_list = pd.read_xml(requests.get(url).text)
      except ValueError: # use all proteins associated with the taxon
        return pd.DataFrame()

    # remove the namespace from the columns
    proteome_list.columns = [x.replace('{http://uniprot.org/proteome}', '') for x in proteome_list.columns]
    return proteome_list

  def get_next_link(self, headers):
    """
    UniProt will provide a link to the next batch of proteins.
    We can use a regular expression to extract the URL from the header.

    Args:
      headers (dict): Headers from UniProt API response.
    """
    re_next_link = re.compile(r'<(.+)>; rel="next"') # regex to extract URL
    if 'Link' in headers:
      match = re_next_link.match(headers['Link'])
      if match:
        return match.group(1)

  def get_protein_batches(self, batch_url):
    """
    Get a batches of proteins from UniProt API because it limits the
    number of proteins you can get at once. Yield each batch until the 
    URL link is empty.
    
    Args:
      batch_url (str): URL to get all proteins for a species.
    """
    while batch_url:
      response = requests.get(batch_url)
      response.raise_for_status()
      yield response
      batch_url = self.get_next_link(response.headers)

  def get_all_proteins(self):
    """
    Get every protein associated with a taxon ID on UniProt.
    Species on UniProt will have a proteome, but not every protein is
    stored within those proteomes. There is a way to get every protein
    using the taxonomy part of UniProt. 
    """
    # URL link to all proteins for a species - size = 500 proteins at a time
    url = f'https://rest.uniprot.org/uniprotkb/search?format=fasta&'\
          f'query=taxonomy_id:{self.taxon_id}&size=500' 

    # loop through all protein batches and write proteins to FASTA file
    for batch in self.get_protein_batches(url):
      with open(f'species/{self.taxon_id}/proteome.fasta', 'a') as f:
        f.write(batch.text)

  def get_gp_proteome_to_fasta(self, proteome_id, proteome_taxon):
    """
    Write the gene priority proteome to a file. 
    This is only for representative and reference proteomes.

    Args:
      proteome_id (str): Proteome ID.
      proteome_taxon (str): Taxon ID for the proteome.
    """
    group = self.species_df[self.species_df['Taxon ID'].astype(str) == self.taxon_id]['Group'].iloc[0]
    ftp_url = f'https://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/reference_proteomes/'
    
    if group == 'archeobacterium':
      ftp_url += f'Archaea/{proteome_id}/{proteome_id}_{proteome_taxon}.fasta.gz'
    elif group == 'bacterium':
      ftp_url += f'Bacteria/{proteome_id}/{proteome_id}_{proteome_taxon}.fasta.gz'
    elif group in ['vertebrate', 'other-eukaryote']:
      ftp_url += f'Eukaryota/{proteome_id}/{proteome_id}_{proteome_taxon}.fasta.gz'
    elif group in ['virus', 'small-virus', 'large-virus']:
      ftp_url += f'Viruses/{proteome_id}/{proteome_id}_{proteome_taxon}.fasta.gz'

    # write the gene priority proteome to a file
    with open(f'species/{self.taxon_id}/gp_proteome.fasta', 'wb') as f:
      f.write(gzip.open(requests.get(ftp_url, stream=True).raw, 'rb').read())

  def get_proteome_to_fasta(self, proteome_id):
    """
    Get the FASTA file for a proteome from UniProt. Either through the
    API or the FTP server. If the proteome is a reference proteome, then
    we also need to get the gene priority proteome, which is a file that 
    contains the best protein record for each gene. We will use this for 
    better gene assignment in the assign_genes.py script.

    If the proteome is not a representative or reference proteome, then
    we will just get the FASTA file from the API.

    Args:
      taxon_id (int):       Taxon ID of species.
      group (str):          Group of species. (e.g. bacterium)
      proteome_id (str):    Proteome ID.
      proteome_taxon (int): Taxon ID of proteome.
      proteome_type (str):  Proteome type. Either: 1. Representative, 
                            2. Reference, 3. Non-redundant, 4. Other
    """
    # create directory for the species if it doesn't exist
    os.makedirs(f'species/{self.taxon_id}', exist_ok=True)

    url = f'https://rest.uniprot.org/uniprotkb/stream?query=proteome:{proteome_id}&format=fasta&compressed=false&includeIsoform=true'
    with open(f'species/{self.taxon_id}/{proteome_id}.fasta', 'w') as f:
      f.write(requests.get(url).text)

  def get_proteome_with_most_matches(self, epitopes_df):
    """
    Get the proteome ID and true taxon associated with
    the proteome with the most proteins in case there is a tie.
    We get the true taxon so we can extract the data from the FTP server.
    """
    # if there is only one proteome, then we don't need to do anything else
    # just get the proteome and return the ID and taxon
    if self.num_of_proteomes == 1:
      self.get_proteome_to_fasta(self.proteome_list['upid'].iloc[0])
      return self.proteome_list['upid'].iloc[0], self.proteome_list['taxonomy'].iloc[0]

    # TODO: be able to check for discontinuous epitopes
    epitopes = [epitope for epitope in epitopes_df['Peptide'] if not any(char.isdigit() for char in epitope)]
    
    match_counts = {}
    for proteome_id in list(self.proteome_list['upid']):
      self.get_proteome_to_fasta(proteome_id)
      Preprocessor(f'species/{self.taxon_id}/{proteome_id}.fasta', 'sql', f'species/{self.taxon_id}/').preprocess(k=5)
      matches_df = Matcher(epitopes, proteome_id, 0, 5, f'species/{self.taxon_id}/', output_format='dataframe').match()
      
      matches_df.drop_duplicates(subset=['Query Sequence'], inplace=True)
      match_counts[proteome_id] = matches_df['Matched Sequence'].dropna().count()

    proteome_id = max(match_counts, key=match_counts.get)
    proteome_taxon = self.proteome_list[self.proteome_list['upid'] == proteome_id]['taxonomy'].iloc[0]

    return proteome_id, proteome_taxon

  def remove_other_proteomes(self, proteome_id):
    """
    Remove proteome .fasta and .db files that are not the chosen proteome for that species.
    Also, rename the chosen proteome to proteome.fasta and proteome.db.
    """
    proteome_list_to_remove = self.proteome_list[self.proteome_list['upid'] != proteome_id]
    for i in list(proteome_list_to_remove['upid']):
      os.remove(f'species/{self.taxon_id}/{i}.fasta')
      os.remove(f'species/{self.taxon_id}/{i}.db')
    
    # rename the chosen proteome to proteome.fasta and remove the .db file
    os.rename(f'species/{self.taxon_id}/{proteome_id}.fasta', f'species/{self.taxon_id}/proteome.fasta')
    if self.num_of_proteomes > 1: # there is only a .db file if there is more than one proteome
      os.remove(f'species/{self.taxon_id}/{proteome_id}.db')

  def select_proteome(self, epitopes_df):
    """
    Select the proteome to use for a species and its taxon ID.
    
    Check UniProt proteome file for all proteomes associated with that
    taxon. Then do the following checks:

    1. Are there any representative proteomes?
    2. Are there any reference proteomes?
    3. Are there any non-redudant proteomes?
    4. Are there any other proteomes?

    If no to all of the above, then get every protein associated with
    the taxon using the get_all_proteins method from uniprot.org/taxonomy.

    If yes to any of the above, get the proteome ID and type. 
    """
    # if there is no proteome_list, get all proteins associated with that taxon ID
    if self.proteome_list.empty:
      self.get_all_proteins(self.taxon_id)
      return 'None', self.taxon_id, 'All-proteins'

    if self.proteome_list['isRepresentativeProteome'].any():
      proteome_type = 'Representative'
      self.proteome_list = self.proteome_list[self.proteome_list['isRepresentativeProteome']]
      proteome_id, proteome_taxon = self.get_proteome_with_most_matches(epitopes_df)
      self.get_gp_proteome_to_fasta(proteome_id, proteome_taxon)
    
    elif self.proteome_list['isReferenceProteome'].any():
      proteome_type = 'Reference'
      self.proteome_list = self.proteome_list[self.proteome_list['isReferenceProteome']]
      proteome_id, proteome_taxon = self.get_proteome_with_most_matches(epitopes_df)
      self.get_gp_proteome_to_fasta(proteome_id, proteome_taxon)

    elif 'redundantTo' not in self.proteome_list.columns:
      proteome_type = 'Other'
      proteome_id, proteome_taxon = self.get_proteome_with_most_matches(epitopes_df)
    
    elif self.proteome_list['redundantTo'].isna().any():
      proteome_type = 'Non-redundant'
      self.proteome_list = self.proteome_list[self.proteome_list['redundantTo'].isna()]
      proteome_id, proteome_taxon = self.get_proteome_with_most_matches(epitopes_df)
    
    else:
      proteome_type = 'Other'
      proteome_id, proteome_taxon = self.get_proteome_with_most_matches(epitopes_df)
    
    self.remove_other_proteomes(proteome_id)

    return proteome_id, proteome_taxon, proteome_type

  def proteome_to_csv(self):
    """
    Write the proteome data for a species to a CSV file for later use.
    """
    # read in the FASTA file and then get the gene priority IDs if they exist
    proteins = list(SeqIO.parse(f'species/{self.taxon_id}/proteome.fasta', 'fasta'))
    if os.path.isfile(f'species/{self.taxon_id}/gp_proteome.fasta'):
      gp_ids = [str(protein.id.split('|')[1]) for protein in list(SeqIO.parse('9606_gp.fasta', 'fasta'))]

    # start collecting proteome data
    proteome_data = []
    for protein in proteins:
      uniprot_id = protein.id.split('|')[1]
      gp = 1 if uniprot_id in gp_ids else 0

      # TODO: look into using HUGO to map old gene names to new ones

      try:
        gene = re.search('GN=(.*?) ', protein.description).group(1)
      except AttributeError:
        try: # some gene names are at the end of the header
          gene = re.search('GN=(.*?)$', protein.description).group(1)
        except AttributeError:
          gene = ''
      try: # protein existence level
        pe_level = int(re.search('PE=(.*?) ', protein.description).group(1))
      except AttributeError:
        pe_level = 0
      
      proteome_data.append([protein.id.split('|')[0], gene, uniprot_id, gp, pe_level, str(protein.seq)])
      
    pd.DataFrame(data, columns=['db', 'gene', 'id', 'gp', 'pe_level', 'seq']).to_csv(f'species/{self.taxon_id}/proteome.csv', index=False)

def main():
  # define command line args which will take in a taxon ID, user, and password (for IEDB MySQL connection)
  parser = argparse.ArgumentParser()
  
  parser.add_argument('-u', '--user', required=True, help='User for IEDB MySQL connection.')
  parser.add_argument('-p', '--password', required=True, help='Password for IEDB MySQL connection.')
  parser.add_argument('-t', '--taxon_id', required=True, help='Taxon ID for the species to pull data for.')
  
  args = parser.parse_args()
  user = args.user
  password = args.password
  taxon_id = args.taxon_id

  # read in species file and save all taxon IDs to list for checking
  species_df = pd.read_csv('species.csv')
  valid_taxon_ids = species_df['Taxon ID'].astype(str).tolist()

  # do proteome selection for all IEDB species
  # then call get_proteome_to_fasta for each species
  if taxon_id == 'all':
    proteomes = {}
    for taxon_id in valid_taxon_ids:
      # get data for taxon ID
      Fetcher = DataFetcher(user, password, taxon_id)
      epitopes_df = Fetcher.get_epitopes()
      sources_df = Fetcher.get_sources()

      # select proteome
      Selector = ProteomeSelector(taxon_id)
      num_of_proteomes, proteome_id, proteome_taxon, proteome_type = Selector.select_proteome(epitopes_df)
      # Selector.proteome_to_csv()
      
      proteomes[taxon_id] = (num_of_proteomes, proteome_id, proteome_taxon, proteome_type)
      proteome_to_csv(taxon_id)
      print(f'# of Proteomes: {num_of_proteomes}')
      print(f'Proteome ID: {proteome_id}')
      print(f'Proteome taxon: {proteome_taxon}')
      print(f'Proteome type: {proteome_type}')
    
    # create new species columns with proteome ID and type - save to file
    species_df['# of Proteomes'] = species_df['Taxon ID'].map(proteomes).map(lambda x: x[0])
    species_df['Proteome ID'] = species_df['Taxon ID'].map(proteomes).map(lambda x: x[1])
    species_df['Proteome Taxon'] = species_df['Taxon ID'].map(proteomes).map(lambda x: x[2])
    species_df['Proteome Type'] = species_df['Taxon ID'].map(proteomes).map(lambda x: x[3])

    species_df.to_csv('species.csv', index=False)

  # or just one species at a time - check if its valid
  else:
    assert taxon_id in valid_taxon_ids, f'{taxon_id} is not a valid taxon ID.'
    
    # get data for taxon ID
    Fetcher = DataFetcher(user, password, taxon_id)
    epitopes_df = Fetcher.get_epitopes()
    sources_df = Fetcher.get_sources()

    Selector = ProteomeSelector(taxon_id)
    num_of_proteomes, proteome_id, proteome_taxon, proteome_type = Selector.select_proteome(epitopes_df)
    # Selector.proteome_to_csv()

    print(f'# of Proteomes: {num_of_proteomes}')
    print(f'Proteome ID: {proteome_id}')
    print(f'Proteome taxon: {proteome_taxon}')
    print(f'Proteome type: {proteome_type}')

if __name__ == '__main__':
  main()