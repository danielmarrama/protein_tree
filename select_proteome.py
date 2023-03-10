#!/usr/bin/env python3

import warnings
warnings.filterwarnings('ignore')

import re
import os
import pandas as pd
import requests

from pepmatch import Preprocessor, Matcher


class ProteomeSelector:
  def __init__(self, taxon_id):
    self.species_df = pd.read_csv('species.csv')
    self.taxon_id = taxon_id
    
    # create species path with species taxon and name; example: "24-Shewanella putrefaciens"
    species_id_to_name_map = dict(zip(self.species_df['Taxon ID'].astype(str), self.species_df['Species Label']))
    self.species_path = f'species/{taxon_id}-{species_id_to_name_map[taxon_id].replace(" ", "_")}'

    # get proteome list for species and count number of proteomes
    self.proteome_list = self._get_proteome_list()
    self.num_of_proteomes = len(self.proteome_list) + 1 # +1 because "all proteins" is also a candidate proteome


  def select_proteome(self, epitopes_df):
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
    """
    # if species_dir already exists then return the already selected proteome, else create dir
    if os.path.exists(f'{self.species_path}/proteome.fasta'):
      proteome_id = self.species_df[self.species_df['Taxon ID'].astype(str) == self.taxon_id]['Proteome ID'].iloc[0]
      proteome_type = self.species_df[self.species_df['Taxon ID'].astype(str) == self.taxon_id]['Proteome Type'].iloc[0]
      return proteome_id, self.taxon_id, proteome_type
    else:
      os.makedirs(self.species_path, exist_ok=True)

    # if there is no proteome_list, get all proteins associated with that taxon ID
    if self.proteome_list.empty:
      self._get_all_proteins()
      return 'None', self.taxon_id, 'All-proteins'

    if self.proteome_list['isRepresentativeProteome'].any():
      proteome_type = 'Representative'
      self.proteome_list = self.proteome_list[self.proteome_list['isRepresentativeProteome']]
      proteome_id, proteome_taxon = self._get_proteome_with_most_matches(epitopes_df)
      self._get_gp_proteome_to_fasta(proteome_id, proteome_taxon)
    
    elif self.proteome_list['isReferenceProteome'].any():
      proteome_type = 'Reference'
      self.proteome_list = self.proteome_list[self.proteome_list['isReferenceProteome']]
      proteome_id, proteome_taxon = self._get_proteome_with_most_matches(epitopes_df)
      self._get_gp_proteome_to_fasta(proteome_id, proteome_taxon)

    elif 'redundantTo' not in self.proteome_list.columns:
      proteome_type = 'Other'
      proteome_id, proteome_taxon = self._get_proteome_with_most_matches(epitopes_df)
    
    elif self.proteome_list['redundantTo'].isna().any():
      proteome_type = 'Non-redundant'
      self.proteome_list = self.proteome_list[self.proteome_list['redundantTo'].isna()]
      proteome_id, proteome_taxon = self._get_proteome_with_most_matches(epitopes_df)
    
    else:
      proteome_type = 'Other'
      proteome_id, proteome_taxon = self._get_proteome_with_most_matches(epitopes_df)
    
    self._remove_other_proteomes(proteome_id)

    return proteome_id, proteome_taxon, proteome_type

  def proteome_to_csv(self):
    """
    Write the proteome data for a species to a CSV file for later use.
    """
    from Bio import SeqIO

    # read in the FASTA file and then get the gene priority IDs if they exist
    proteins = list(SeqIO.parse(f'{self.species_path}/proteome.fasta', 'fasta'))

    gp_proteome_path = f'{self.species_path}/gp_proteome.fasta'
    if os.path.isfile(gp_proteome_path):
      gp_ids = [str(protein.id.split('|')[1]) for protein in list(SeqIO.parse(gp_proteome_path, 'fasta'))]
    else:
      gp_ids = []

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
    
    columns = ['Database', 'Gene Symbol', 'UniProt ID', 'Gene Priority', 'Protein Existence Level', 'Sequence']
    pd.DataFrame(proteome_data, columns=columns).to_csv(f'{self.species_path}/proteome.csv', index=False)

  def _get_proteome_list(self):
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

  def _get_all_proteins(self):
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
    for batch in self._get_protein_batches(url):
      with open(f'{self.species_path}/proteome.fasta', 'a') as f:
        f.write(batch.text)

  def _get_protein_batches(self, batch_url):
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

  def _get_next_link(self, headers):
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

  def _get_gp_proteome_to_fasta(self, proteome_id, proteome_taxon):
    """
    Write the gene priority proteome to a file. 
    This is only for representative and reference proteomes.
    Depending on the species group, the FTP URL will be different.

    Args:
      proteome_id (str): Proteome ID.
      proteome_taxon (str): Taxon ID for the proteome.
    """
    import gzip
    
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
    
    r = requests.get(ftp_url, stream=True)
    try:
      r.raise_for_status()
    except:
      return

    # unzip the request and write the gene priority proteome to a file
    with open(f'{self.species_path}/gp_proteome.fasta', 'wb') as f:
      f.write(gzip.open(r.raw, 'rb').read())

  def _get_proteome_to_fasta(self, proteome_id):
    """
    Get the FASTA file for a proteome from UniProt API.
    Include all isoforms and do not compress the file.
    """
    url = f'https://rest.uniprot.org/uniprotkb/stream?compressed=false&format=fasta&includeIsoform=true&query=(proteome:{proteome_id})'
    r = requests.get(url)
    r.raise_for_status()
    with open(f'{self.species_path}/{proteome_id}.fasta', 'w') as f:
      f.write(r.text)

  def _get_proteome_with_most_matches(self, epitopes_df):
    """
    Get the proteome ID and true taxon ID associated with
    that proteome with the most epitope matches in case there is a tie.
    We get the true taxon so we can extract the data from the FTP server
    if needed.
    """
    # if there is only one proteome, then we don't need to do anything else
    # just get the proteome and return the ID and taxon
    if self.num_of_proteomes <= 2:
      self._get_proteome_to_fasta(self.proteome_list['upid'].iloc[0])
      return self.proteome_list['upid'].iloc[0], self.proteome_list['taxonomy'].iloc[0]

    # drop any epitopes that do not have a peptide sequence - this is rare but needs to be checked
    epitopes_df = epitopes_df[epitopes_df['Peptide'].notna()]
    # TODO: be able to check for discontinuous epitopes
    epitopes = [epitope for epitope in epitopes_df['Peptide'] if not any(char.isdigit() for char in epitope)]
    
    # keep track of # of epitope matches for each proteome
    match_counts = {}
    for proteome_id in list(self.proteome_list['upid']):
      self._get_proteome_to_fasta(proteome_id)
      
      # PEPMatch preprocessing and matching
      Preprocessor(f'{self.species_path}/{proteome_id}.fasta', 'sql', f'{self.species_path}').preprocess(k=5)
      matches_df = Matcher(epitopes, proteome_id, 0, 5, f'{self.species_path}', output_format='dataframe').match()
      
      # remove any duplicate matches and count the number of matches
      matches_df.drop_duplicates(subset=['Query Sequence'], inplace=True)
      try: 
        match_counts[proteome_id] = matches_df['Matched Sequence'].dropna().count()
      except KeyError:
        match_counts[proteome_id] = 0 # in case there are no matches

    # select the proteome ID and proteome taxon with the most matches
    proteome_id = max(match_counts, key=match_counts.get)
    proteome_taxon = self.proteome_list[self.proteome_list['upid'] == proteome_id]['taxonomy'].iloc[0]

    return proteome_id, proteome_taxon

  def _remove_other_proteomes(self, proteome_id):
    """
    Remove the proteome FASTA files that are not the chosen proteome for that
    species. Also, remove the .db files and rename the chosen proteome to 
    "proteome.fasta".
    """
    proteome_list_to_remove = self.proteome_list[self.proteome_list['upid'] != proteome_id]
    for i in list(proteome_list_to_remove['upid']):
      os.remove(f'{self.species_path}/{i}.fasta')
      os.remove(f'{self.species_path}/{i}.db')
    
    # rename the chosen proteome to proteome.fasta and remove the .db file
    os.rename(f'{self.species_path}/{proteome_id}.fasta', f'{self.species_path}/proteome.fasta')
    if self.num_of_proteomes > 2: # there is only a .db file if there is more than one proteome
      os.remove(f'{self.species_path}/{proteome_id}.db')


def main():
  import argparse
  from get_data import DataFetcher

  # define command line args which will take in a taxon ID, user, and password (for IEDB MySQL connection)
  parser = argparse.ArgumentParser()
  
  parser.add_argument('-u', '--user', required=True, help='User for IEDB MySQL connection.')
  parser.add_argument('-p', '--password', required=True, help='Password for IEDB MySQL connection.')
  parser.add_argument('-t', '--taxon_id', required=True, help='Taxon ID for the species to pull data for.')
  
  args = parser.parse_args()
  user = args.user
  password = args.password
  taxon_id = args.taxon_id

  # read in IEDB species data
  species_df = pd.read_csv('species.csv')
  valid_taxon_ids = species_df['Taxon ID'].astype(str).tolist()

  # dicts for mapping taxon IDs to all their taxa and their names
  all_taxa_map = dict(zip(species_df['Taxon ID'].astype(str), species_df['All Taxa']))
  species_id_to_name_map = dict(zip(species_df['Taxon ID'].astype(str), species_df['Species Label']))

  # do proteome selection for all IEDB species
  if taxon_id == 'all':
    proteomes = {}
    for taxon_id in valid_taxon_ids:
      # get data for taxon ID
      Fetcher = DataFetcher(user, password, taxon_id, all_taxa_map[taxon_id])
      epitopes_df = Fetcher.get_epitopes()
      sources_df = Fetcher.get_sources()

      # select proteome
      Selector = ProteomeSelector(taxon_id)
      num_of_proteomes, proteome_id, proteome_taxon, proteome_type = Selector.select_proteome(epitopes_df)
      Selector.proteome_to_csv()
      
      proteomes[taxon_id] = (num_of_proteomes, proteome_id, proteome_taxon, proteome_type)
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
    Fetcher = DataFetcher(user, password, taxon_id, all_taxa_map[taxon_id])
    epitopes_df = Fetcher.get_epitopes()
    sources_df = Fetcher.get_sources()

    Selector = ProteomeSelector(taxon_id)
    print(f'Number of candidate proteomes: {Selector.num_of_proteomes}')
    proteome_id, proteome_taxon, proteome_type = Selector.select_proteome(epitopes_df)
    Selector.proteome_to_csv()

    print(f'Proteome ID: {proteome_id}')
    print(f'Proteome taxon: {proteome_taxon}')
    print(f'Proteome type: {proteome_type}')

if __name__ == '__main__':
  main()