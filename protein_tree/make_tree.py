#!/usr/bin/env python3

import math
import pandas as pd
from anytree import Node, RenderTree
from pathlib import Path


def build_tree_for_species(peptide_assignments_df):
  gene_nodes = {}
  protein_tracker = {}

  for _, row in peptide_assignments_df.iterrows():
    gene = row['Assigned Gene'] if not pd.isna(row['Assigned Gene']) else 'Unknown Gene'
    protein_name = row['Assigned Protein Name']
    protein_id = row['Assigned Protein ID']

    if pd.isna(protein_id):
      continue

    protein_info = f"{protein_name} (UniProt: {protein_id})" if not pd.isna(protein_name) else f"UniProt: {protein_id}"

    if gene not in gene_nodes:
      gene_nodes[gene] = Node(gene)
      protein_tracker[gene] = set()

    if protein_info not in protein_tracker[gene]:
      Node(protein_info, parent=gene_nodes[gene])
      protein_tracker[gene].add(protein_info)

  return gene_nodes


if __name__ == '__main__':

  import argparse

  parser = argparse.ArgumentParser()

  parser.add_argument(
    '-b', '--build_path', 
    type=str,
    default=Path(__file__).parent.parent / 'build',
    help='Path to data directory.'
  )
  parser.add_argument(
    '-a', '--all_species', 
    action='store_true', 
    help='Build protein tree for all IEDB species.'
  )
  parser.add_argument(
    '-t', '--taxon_id', 
    type=int,
    help='Taxon ID for the species to pull data for.'
  )

  args = parser.parse_args()

  build_path = Path(args.build_path)
  all_species = args.all_species
  taxon_id = args.taxon_id

  if all_species:
    pass

  else:
    peptide_assignments_df = pd.read_csv(build_path / 'species' / str(taxon_id) / 'peptide-assignments.tsv', sep='\t')
    peptide_assignments_df['Assigned Gene'].fillna(peptide_assignments_df['ARC Assignment'], inplace=True)
    tree = build_tree_for_species(peptide_assignments_df)

    for nodes in tree:
      for pre, _, node in RenderTree(tree[nodes]):
        print("%s%s" % (pre, node.name))
