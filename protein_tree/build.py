#!/usr/bin/env python3

import pandas as pd
import sqlite3
from argparse import ArgumentParser



def build_old_tree(df):
  pass


def build_new_tree(df):
  pass




def main():
  parser = ArgumentParser('Build the protein tree.')
 
  parser.add_argument('sqlite_db', help='Path to the SQLite database')
  # parser.add_argument(
  #   'peptide_assignments',
  #   help='Path to peptide-assignments.tsv to read',
  # )
  args = parser.parse_args()

  with sqlite3.connect(args.sqlite_db) as connection:
    df = pd.read_sql_query("SELECT * FROM organism_tree", connection)
  
  predicates = ['rdfs:subClassOf', 'rdf:type', 'rdfs:label']
  lower_subjects = df[(df['predicate'] == 'iedb-taxon:level') & (df['object'] == 'lower')]['subject']

  df = df[(df['predicate'].isin(predicates))]
  df = df[~(df['subject'].isin(lower_subjects))]
  df['graph'] = 'iedb-taxon:protein_tree'
  df.loc[df['predicate'] == 'rdfs:label', 'object'] = df['object'] + ' protein'
  df.loc[df['object'] == 'Organism protein', 'object'] = 'protein'
  
  build_old_tree(df)
  build_new_tree(df)
  df.to_csv('~/Downloads/protein_tree.csv', index=False)

if __name__ == "__main__":
  main()
