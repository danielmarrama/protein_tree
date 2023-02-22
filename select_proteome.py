#!/usr/bin/env python3

import re
import requests

from requests.adapters import HTTPAdapter, Retry

from Bio import SeqIO


# class ProteomeSelector(object):
#   def __init__(self, taxon_id):
#     self.taxon_id = taxon_id



def get_next_link(headers):
  """
  UniProt will provide a link to the next batch of proteins.
  We can use a regular expression to extract the URL from the header.

  Args:
    headers (dict): Headers from UniProt API response.
  """
  re_next_link = re.compile(r'<(.+)>; rel="next"')
  if "Link" in headers:
    match = re_next_link.match(headers["Link"])
    if match:
      return match.group(1)

def get_protein_batches(batch_url):
  """
  Get a batches of proteins from UniProt API because it limits the
  number of proteins you can get at once. Yield each batch until the 
  next link is empty.
  
  Args:
    batch_url (str): URL to get.
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

  Args:
    taxon_id (str): Taxon ID to select.
  """
  # URL link to all proteins for a species - size = 500 proteins at a time
  url = f'https://rest.uniprot.org/uniprotkb/search?format=fasta&'\
        f'query=taxonomy_id:{taxon_id}&size=500' 

  for batch in get_protein_batches(url):
    with open(f'{taxon_id}.fasta', 'a') as f:
      f.write(batch.text)

if __name__ == '__main__':
  get_all_proteins('9606')