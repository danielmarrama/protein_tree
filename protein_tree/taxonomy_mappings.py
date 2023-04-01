#!/usr/bin/env python3

import requests
import _pickle as pickle


def get_taxonomy_information(taxon_id):
  url = f"https://rest.uniprot.org/taxonomy/{taxon_id}.json"
  response = requests.get(url)
  taxonomy_data = response.json()
  return taxonomy_data

def get_species_data(taxon_id, taxonomy_data):
  try:
    if taxonomy_data["rank"] in ["species", "species group", "genus", "family"]:
      return taxonomy_data["taxonId"], taxonomy_data["scientificName"]
    else:
      for lineage_item in taxonomy_data["lineage"]:
        if lineage_item["rank"] == "species":
          return lineage_item["taxonId"], lineage_item["scientificName"]
  except:
      return taxon_id, ''

def create_species_mapping(taxon_ids):
  species_mapping = {}
  for taxon_id in taxon_ids:
    print(taxon_id)
    taxonomy_data = get_taxonomy_information(taxon_id)
    species_id, species_name = get_species_data(taxon_id, taxonomy_data)
    print(species_id, species_name)
    species_mapping[taxon_id] = (species_id, species_name)

  # os.path.realpath(__file__)
  with open('species_mapping.pickle', 'wb') as f:
    pickle.dump(species_mapping, f)

  return species_mapping


# def get_relevant_organisms(user, password) -> pd.DataFrame:
#   sql_query = """
#               SELECT *
#               FROM organism
#               WHERE RANK in ("species", "subspecies", "varietas", "forma", "no rank");
#               """
#   sql_engine = create_sql_engine(user, password)
#   return pd.DataFrame(sql_engine.connect().execute(text(sql_query)))

# def lower_rank_to_species_map(df, tax_id_to_row, lower_ranks) -> dict:
#   """
#   Create a dictionary mapping of all lower ranks to its parent species rank.
#   """
#   lower_rank_to_species_mapping = {}
#   for _, row in df.iterrows():
#     if row['rank'] in lower_ranks:
#       # get the path and split it into a list of taxon IDs
#       path_list = row['path'].split(':')

#       # iterate through the parent_tax_id_string in reverse order (from lowest to highest rank)
#       for parent_tax_id in reversed(path_list):

#         # find the row in the dictionary corresponding to the current taxon ID
#         parent_row = tax_id_to_row.get(int(parent_tax_id))
#         if parent_row is not None:
#           # check if the current taxon ID is a species
#           if parent_row['rank'] == 'species':
#             # map the current row's taxon ID to the species taxon ID
#             lower_rank_to_species_mapping[str(row['tax_id'])] = (str(parent_row['tax_id']), str(parent_row['organism_name']))
#             break
#   return lower_rank_to_species_mapping

# def species_to_lower_ranks_map(df, tax_id_to_row, lower_ranks) -> dict:
#   """
#   Create a dictionary mapping of species rank to all lower ranks.
#   """
#   species_to_lower_ranks_mapping = {}
#   for _, row in df.iterrows():
#     if row['rank'] in lower_ranks:
      
#       # get the taxonomy path and split it into a list of taxon IDs
#       path_list = row['path'].split(':')

#       # iterate through the parent_tax_id_string in reverse order (from lowest to highest rank)
#       for parent_tax_id in reversed(path_list):

#         # check if the parent taxon is a species
#         parent_row = tax_id_to_row.get(int(parent_tax_id))
#         if parent_row is not None and parent_row['rank'] == 'species':
           
#             # create a new key in the dictionary if it doesn't exist
#             if str(parent_row['tax_id']) not in species_to_lower_ranks_mapping:
#               species_to_lower_ranks_mapping[str(parent_row['tax_id'])] = []
            
#             # check if the taxon ID is not 'None' before appending it to the list
#             if row["tax_id"] is not None and str(row["tax_id"]) != 'None':
#                 species_to_lower_ranks_mapping[str(parent_row["tax_id"])].append(str(row["tax_id"]))
#             break
#   return species_to_lower_ranks_mapping

# def species_id_to_name_map(df) -> dict:
#   """
#   Create a dictionary mapping of species taxon ID to its name.
#   """
#   df = df[df['rank'] == 'species']
#   return dict(zip(df['tax_id'].astype(str), df['organism_name']))

# def create_mappings(user, password):
#   """
#   Create dicts of mappings from lower ranks to species and vice versa.
#   Store them in pickle files.
#   """
#   # create mappings directory if it doesn't exist
#   os.makedirs('mappings', exist_ok=True)

#   # retrieve organism table with species and all other lower ranks
#   df = get_relevant_organisms(user, password)

#   tax_id_to_row = {row["tax_id"]: row for _, row in df.iterrows()}
#   lower_ranks = ["subspecies", "varietas", "forma", "no rank"]

#   lower_rank_to_species_mapping  = lower_rank_to_species_map(df, tax_id_to_row, lower_ranks)
#   species_to_lower_ranks_mapping = species_to_lower_ranks_map(df, tax_id_to_row, lower_ranks)
#   species_id_to_name_mapping     = species_id_to_name_map(df)

#   with open('mappings/lower_rank_to_species_mapping.pickle', 'wb') as f:
#     pickle.dump(lower_rank_to_species_mapping, f)

#   with open('mappings/species_to_lower_ranks_mapping.pickle', 'wb') as f:
#     pickle.dump(species_to_lower_ranks_mapping, f)

#   with open('mappings/species_id_to_name_mapping.pickle', 'wb') as f:
#     pickle.dump(species_id_to_name_mapping, f)


# if __name__ == '__main__':
#   import argparse
#   parser = argparse.ArgumentParser()
#   parser.add_argument('-u', '--user', required=True, help='MySQL backend username.')
#   parser.add_argument('-p', '--password', required=True, help='MySQL backend password.')
#   args = parser.parse_args()

#   create_mappings(args.user, args.password)