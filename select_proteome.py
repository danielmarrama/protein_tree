#!/usr/bin/env python3

import warnings
warnings.filterwarnings('ignore')

import io
import re
import os
import argparse
import pandas as pd
import numpy as np

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
    with open(f'species/{taxon_id}/{taxon_id}.fasta', 'a') as f:
      f.write(batch.text)

def get_proteome_to_fasta(taxon_id, proteome_id, proteome_type):
  """
  Get the FASTA file for a proteome from UniProt. Either through the
  API or the FTP server. If the proteome is a reference proteome, then
  we also need to get the gene priority proteome, which is a file that 
  contains the best protein record for each gene. We will use this for 
  better gene assignment in the assign_genes.py script.

  Args:
    proteome_id (str):   Proteome ID.
    proteome_type (str): Proteome type. Either: 1. Representative, 
                         2. Reference, 3. Non-redundant, 4. Other
  """
  # the species without a formal proteome need to be skipped
  # the proteins for them are extracted using get_all_proteins
  if proteome_type == 'All-proteins':
    get_all_proteins(taxon_id)
    return

  os.makedirs(f'species/{taxon_id}', exist_ok=True)

  url = f'https://rest.uniprot.org/uniprotkb/stream?query=proteome:{proteome_id}&format=fasta&compressed=false&includeIsoform=true'
  # with open()

def get_id_with_max_proteins(proteome_list):
  """Get the proteome ID with the most proteins in case there is a tie."""
  return proteome_list.loc[proteome_list['proteinCount'].idxmax()]['upid']

def select_proteome(taxon_id):
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
  # URL to get proteome list for a species - use proteome_type:1 first to get reference proteomes
  url = f'https://rest.uniprot.org/proteomes/stream?format=xml&query=(proteome_type:1)AND(taxonomy_id:{taxon_id})'
  
  # read in proteome list for the species from the UniProt API - try to get
  # the proteome list with proteome_type:1, which are the reference proteomes
  try:
    proteome_list = pd.read_xml(requests.get(url).text)
  except ValueError:
    try:
      url = url.replace('(proteome_type:1)AND', '')
      proteome_list = pd.read_xml(requests.get(url).text)
    except ValueError:
      return '', 'All-proteins'

  # remove the namespace from the columns
  proteome_list.columns = [x.replace('{http://uniprot.org/proteome}', '') for x in proteome_list.columns]

  # get proteome ID and proteome type based on the checks - if there are ties
  # then get the ID with most proteins
  if proteome_list['isRepresentativeProteome'].any():
    proteome_list = proteome_list[proteome_list['isRepresentativeProteome']]
    return get_id_with_max_proteins(proteome_list), 'Representative'
  elif proteome_list['isReferenceProteome'].any():
    proteome_list = proteome_list[proteome_list['isReferenceProteome']]
    return get_id_with_max_proteins(proteome_list), 'Reference'
  
  if 'redundantTo' not in proteome_list.columns:
    return get_id_with_max_proteins(proteome_list), 'Other'
  elif proteome_list['redundantTo'].isna().any():
    proteome_list = proteome_list[proteome_list['redundantTo'].isna()]
    return get_id_with_max_proteins(proteome_list), 'Non-redundant'
  else:
    return get_id_with_max_proteins(proteome_list), 'Other'

def main():
  # define command line args which will take in a taxon ID
  parser = argparse.ArgumentParser(description='Select and get the proteome to use for a species.')
  parser.add_argument('taxon_id', help='Taxon ID for the species.')
  args = parser.parse_args()
  taxon_id = args.taxon_id

  # read in species file and save all taxon IDs to list for checking
  species_df = pd.read_csv('species.csv')
  valid_taxon_ids = species_df['Taxon ID'].astype(str).tolist()

  # make species folder to save proteome FASTA files
  os.makedirs('species', exist_ok=True)

  # do proteome selection for all IEDB species
  # then call get_proteome_to_fasta for each species
  if taxon_id == 'all':
    proteomes = {}
    for taxon_id in species_df['Taxon ID']:
      proteome, proteome_type = select_proteome(taxon_id)
      proteomes[taxon_id] = (proteome, proteome_type)
      # get_proteome_to_fasta(taxon_id, proteome, proteome_type)
      print(taxon_id, proteome, proteome_type)
    
    # create new species columns with proteome ID and type - save to file
    species_df['Proteome ID'] = species_df['Taxon ID'].map(proteomes).map(lambda x: x[0])
    species_df['Proteome Type'] = species_df['Taxon ID'].map(proteomes).map(lambda x: x[1])
    species_df.to_csv('species.csv', index=False)

  # or just one species at a time - check if its valid
  else:
    assert taxon_id in valid_taxon_ids, f'{taxon_id} is not a valid taxon ID.'
    proteome, proteome_type = select_proteome(taxon_id)
    # get_proteome_to_fasta(taxon_id, proteome, proteome_type)
    print(taxon_id, proteome, proteome_type)

if __name__ == '__main__':
  main()