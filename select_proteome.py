#!/usr/bin/env python3

import warnings
warnings.filterwarnings('ignore')

import re
import argparse
import pandas as pd

import requests
from requests.adapters import HTTPAdapter, Retry

from Bio import SeqIO
from pepmatch import Preprocessor, Matcher


def get_next_link(headers):
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

def get_protein_batches(batch_url):
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
    batch_url = get_next_link(response.headers)

def get_all_proteins(taxon_id):
  """
  Get every protein associated with a taxon ID on UniProt.
  Species on UniProt will have a proteome, but not every protein is
  stored within those proteomes. There is a way to get every protein
  using the taxonomy part of UniProt. 
  """
  # URL link to all proteins for a species - size = 500 proteins at a time
  url = f'https://rest.uniprot.org/uniprotkb/search?format=fasta&'\
        f'query=taxonomy_id:{taxon_id}&size=500' 

  # loop through all protein batches and write proteins to FASTA file
  for batch in get_protein_batches(url):
    with open(f'{taxon_id}.fasta', 'a') as f:
      f.write(batch.text)

def select_proteome(taxon_id):
  """
  Select the proteome to use for a species and its taxon ID.
  
  Check UniProt proteome file for all proteomes associated with that
  taxon. Then do the following checks:

  1. Are there any representative proteomes?
  2. Are there any reference proteomes?
  3. Are there any non-redudant proteome?

  If no to all of the above, then get every protein associated with
  the taxon using the get_all_proteins method from uniprot.org/taxonomy.

  If yes to any of the above, get the proteome ID and download from the
  UniProt FTP server.
  """
  # get list of  proteomes for the species using the taxon ID
  taxon_df = proteomes[proteomes['speciesTaxon'] == int(taxon_id)]

  # make sure there are proteomes for the species
  if taxon_df.empty:
    print(f'No proteomes for {taxon_id}.')
    # get_all_proteins(taxon_id)

  # check if there is a representative proteome
  if taxon_df['isRepresentativeProteome'].any():
    proteome_id = taxon_df[taxon_df['isRepresentativeProteome']]['upid'].iloc[0]
  elif taxon_df['isReferenceProteome'].any():
    proteome_id = taxon_df[taxon_df['isReferenceProteome']]['upid'].iloc[0]
  else:
    proteome_id = str(len(taxon_df))

  return proteome_id

def main():
  # read in proteomes file
  global proteomes 
  proteomes = pd.read_csv('proteomes.csv')
  
  # define command line args which will take in a taxon ID
  parser = argparse.ArgumentParser(description='Select and get the proteome to use for a species.')
  parser.add_argument('taxon_id', help='Taxon ID for the species.')
  args = parser.parse_args()

  taxon_id = args.taxon_id

  # TODO: add option select the proteomes for all species
  if taxon_id == 'all':
    species_df = pd.read_csv('species.csv')
    for taxon_id in species_df['Taxon ID']:
      proteome_id = select_proteome(taxon_id)
      print(taxon_id, proteome_id)

  # proteome_id = select_proteome(taxon_id)

if __name__ == '__main__':
  main()