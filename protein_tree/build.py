#!/usr/bin/env python3

import pandas as pd
import sqlite3
from argparse import ArgumentParser

import warnings
warnings.filterwarnings('ignore')



def build_old_tree(tree_df, peptide_assignments):
  new_rows = []
  for _, row in peptide_assignments.iterrows():
    new_rows.append(old_protein_id(row))
    new_rows.append(old_protein_label(row))
    new_rows.append(old_gene_label(row))
  
  tree_df = pd.concat([tree_df, pd.DataFrame(new_rows)], ignore_index=True)
  
  return tree_df

def old_protein_id(row):
  protein_label_row = {
    'assertion': 1,
    'retraction': 0,
    'graph': 'iedb-taxon:protein_tree',
    'subject': f"{row['Assigned Protein ID']}",
    'predicate': 'rdfs:subClassOf',
    'object': f"NCBITaxon:{row['Species Taxon ID']}",
    'datatype': 'xsd:string',
    'annotation': None
  }
  return protein_label_row

def old_protein_label(row):
  protein_label_row = {
    'assertion': 1,
    'retraction': 0,
    'graph': 'iedb-taxon:protein_tree',
    'subject': f"{row['Assigned Protein ID']}",
    'predicate': 'rdfs:label',
    'object': f"{row['Assigned Protein Name']} (UniProt:{row['Assigned Protein ID']})",
    'datatype': 'xsd:string',
    'annotation': None
  }
  return protein_label_row

def old_gene_label(row):
  gene_label_row = {
    'assertion': 1,
    'retraction': 0,
    'graph': 'iedb-taxon:protein_tree',
    'subject': f"{row['Assigned Protein ID']}",
    'predicate': 'from_gene',
    'object': f"{row['Assigned Gene']}",
    'datatype': 'xsd:string',
    'annotation': None
  }
  return gene_label_row

def build_new_tree(tree_df, peptide_assignments):
  new_rows = []
  for _, row in peptide_assignments.iterrows():
    new_rows.extend(new_gene_label(row))
    new_rows.append(new_protein_id(row))
    new_rows.append(new_protein_label(row))
    
  tree_df = pd.concat([tree_df, pd.DataFrame(new_rows)], ignore_index=True)
  
  return tree_df

def new_gene_label(row):
  gene_class_row = {
    'assertion': 1,
    'retraction': 0,
    'graph': 'iedb-taxon:protein_tree',
    'subject': f"{row['Species Taxon ID']}:{row['Assigned Gene']}",
    'predicate': 'rdf:type',
    'object': 'owl:Class',
    'datatype': '_IRI',
    'annotation': None
  }
  gene_label_row = {
    'assertion': 1,
    'retraction': 0,
    'graph': 'iedb-taxon:protein_tree',
    'subject': f"{row['Species Taxon ID']}:{row['Assigned Gene']}",
    'predicate': 'rdfs:label',
    'object': f"{row['Assigned Gene']}",
    'datatype': 'xsd:string',
    'annotation': None
  }
  gene_subclass_row = {
    'assertion': 1,
    'retraction': 0,
    'graph': 'iedb-taxon:protein_tree',
    'subject': f"{row['Species Taxon ID']}:{row['Assigned Gene']}",
    'predicate': 'rdfs:subClassOf',
    'object': f"NCBITaxon:{row['Species Taxon ID']}",
    'datatype': 'xsd:string',
    'annotation': None
  }
  return gene_class_row, gene_label_row, gene_subclass_row

def new_protein_id(row):
  protein_label_row = {
    'assertion': 1,
    'retraction': 0,
    'graph': 'iedb-taxon:protein_tree',
    'subject': f"{row['Assigned Protein ID']}",
    'predicate': 'rdfs:subClassOf',
    'object': f"{row['Species Taxon ID']}:{row['Assigned Gene']}",
    'datatype': 'xsd:string',
    'annotation': None
  }
  return protein_label_row

def new_protein_label(row):
  protein_label_row = {
    'assertion': 1,
    'retraction': 0,
    'graph': 'iedb-taxon:protein_tree',
    'subject': f"{row['Assigned Protein ID']}",
    'predicate': 'rdfs:label',
    'object': f"{row['Assigned Protein Name']} (UniProt:{row['Assigned Protein ID']})",
    'datatype': 'xsd:string',
    'annotation': None
  }
  return protein_label_row

def main():
  parser = ArgumentParser('Build the protein tree.')
 
  parser.add_argument('sqlite_db', help='Path to the SQLite database')
  parser.add_argument(
    'peptide_assignments',
    help='Path to all-peptide-assignments.tsv to read for tree building.',
  )
    
  args = parser.parse_args()

  peptide_assignments = pd.read_csv(args.peptide_assignments, sep='\t')
  peptide_assignments['Assigned Gene'].fillna(peptide_assignments['ARC Assignment'], inplace=True)
  peptide_assignments.drop_duplicates(subset=['Assigned Protein ID'], inplace=True)

  with sqlite3.connect(args.sqlite_db) as connection:
    tree_df = pd.read_sql_query("SELECT * FROM organism_tree", connection)
  
    predicates = ['rdfs:subClassOf', 'rdf:type', 'rdfs:label']
    lower_subjects = tree_df[tree_df['object'].isin(['upper', 'species'])]['subject']

    tree_df = tree_df[(tree_df['predicate'].isin(predicates))]
    tree_df = tree_df[tree_df['subject'].isin(lower_subjects)]
    tree_df['graph'] = 'iedb-taxon:protein_tree'
    tree_df.loc[tree_df['predicate'] == 'rdfs:label', 'object'] = tree_df['object'] + ' protein'
    tree_df.loc[tree_df['object'] == 'Organism protein', 'object'] = 'protein'
    
    old_df = build_old_tree(tree_df, peptide_assignments)
    new_df = build_new_tree(tree_df, peptide_assignments)
  
    old_df.to_sql('protein_tree_old', connection, if_exists='replace', index=False)
    new_df.to_sql('protein_tree_new', connection, if_exists='replace', index=False)

if __name__ == "__main__":
  main()
