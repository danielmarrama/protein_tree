#!/usr/bin/env python3

import re
import pickle
import pandas as pd

from Bio import SeqIO
from anytree import Node, RenderTree


def pickle_gp_mapping(df):
  '''
  Creates a map of the canonical UniProt protein ID for a gene to all the
  other UniProt protein IDs. Dictionary is stored as a pickle file.
  '''
  gp_mapping = {}
  for i, group in df.groupby('gene'):
    for id in list(group['id']):
      try:
        gp_mapping[id] = str(group[group['gp'] == 1]['id'].iloc[0])
      except IndexError:
        pass

  with open('gp_mapping.pickle', 'wb') as f:
    pickle.dump(gp_mapping, f)


def add_nodes(nodes, parent, child):
  '''Creates nodes for gene to UniProt ID mappings.'''
  if parent not in nodes:
    nodes[parent] = Node(parent)  
  if child not in nodes:
    nodes[child] = Node(child)
  nodes[child].parent = nodes[parent]


def create_protein_tree(proteome):
  '''
  Makes the basic protein tree. Reads in a proteome and an accompanying
  gene priority proteome. Genes as the root and UniProt IDs as the children.
  Outputs a dataframe of the proteome data as well as a text file of the
  protein tree itself. It will also output a pickle file mapping of the 
  canonical UniProt protein ID for a gene to the other UniProt protein IDs.  
  '''
  proteins = list(SeqIO.parse(proteome, 'fasta'))

  # get UniProt IDs for one protein per gene proteome
  gp_ids = [str(x.id.split('|')[1]) for x in list(SeqIO.parse(gp_proteome, 'fasta'))]
  
  data = []
  for protein in proteins:
    uniprot_id = protein.id.split('|')[1]

    # get gene symbol from FASTA file
    try:
      gene = re.search('GN=(.*?) ', protein.description).group(1)
    except AttributeError:

      # sometimes the gene symbol is at the end of the FASTA description
      try:
        gene = re.search('GN=(.*?)$', protein.description).group(1)
      except AttributeError:
        gene = ''

    gp = 1 if uniprot_id in gp_ids else 0
    data.append([protein.id.split('|')[0], gene, uniprot_id, gp, str(protein.seq)])
  
  # put protein tree data into dataframe
  df = pd.DataFrame(data, columns=['db', 'gene', 'id', 'gp', 'seq'])

  # start tree with nodes - genes as root and UniProt IDs as children
  nodes = {}
  for parent, child in zip(df['gene'],df['id']):
    add_nodes(nodes, parent, child)

  # write the tree into a text file
  with open('protein_tree.txt', 'w') as f:
    roots = list(df[~df['gene'].isin(df['id'])]['gene'].unique())
    for root in roots:         # you can skip this for roots[0], if there is no forest and just 1 tree
      for pre, _, node in RenderTree(nodes[root]):
        if node.name in gp_ids:
          f.write("%s%s*" % (pre, node.name))
          f.write('\n')
        else:
          f.write("%s%s" % (pre, node.name))
          f.write('\n')

  pickle_gp_mapping(df)
  df.to_csv('human_proteome.csv', index=False)

  return df