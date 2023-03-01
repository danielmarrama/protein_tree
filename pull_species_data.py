import mysql.connector
import argparse


def get_epitopes(taxon_id):
  """
  Get all epitopes for a species from the IEDB API.
  """
  pass
  # TODO: pull data from the MySQL database


def get_sources(taxon_id):
  """
  Get all source antigens for a species from the IEDB API.
  """
  pass


def main():
  # define command line args which will take in a taxon ID
  parser = argparse.ArgumentParser()
  
  parser.add_argument('-u', '--user', help='User for IEDB MySQL connection.')
  parser.add_argument('-p', '--password', help='User for IEDB MySQL connection.')
  
  args = parser.parse_args()
  user = args.user
  password = args.password

  connection = mysql.connector.connect(
    host="iedb-mysql.liai.org",
    port=33306,
    user=user,
    passwd=password,
    database='iedb_query'
  )

  c = connection.cursor()
  c.execute("SHOW TABLES;")
  print(c.fetchall())

if __name__ == '__main__':
  main()
