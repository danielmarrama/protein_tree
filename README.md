# IEDB Protein Tree 

Mapping epitopes and source antigens to their parent protein (and gene).

### Process
1. Fetch epitopes and source antigens for species in the IEDB.
2. Select the best proteome for that species.
3. Assign genes to source antigens using BLAST and PEPMatch.
4. Assign epitopes to their parent protein.

### Inputs
- List of IEDB species: [species.csv](species.csv)
- List of manual parents: [manual_parents.csv](manual_parents.csv)
    - This is the list of sources that have been manually assigned for their parents.
- IEDB MySQL Backend username and password
- `blastp` and `makeblastdb` executables from [NCBI](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
- Taxon ID for species to build or `-a`(all) to build all species.

### Running
For one species:
``` bash
./protein_tree.py -u <username> -p <password> -t <taxon ID>
```
or for all species:
``` bash
./protein_tree.py -u <username> -p <password> -a
```

### Outputs

For each species:
- proteome.fasta - selected proteome in FASTA
- [optional] gp_proteome.fasta - the gene priority proteome in FASTA if it exists
- proteome.csv   - selected proteome with data in CSV
- sources.csv    - each source antigen with assigned gene
- epitopes.csv   - each epitope with its source antigen and assigned parent protein

For all species:
- metrics.csv    - the metadata: # of source antigens, # of epitopes, % of sources assigned a gene, % of epitopes assigned a parent, etc.

Use [combine_data.py](combine_data.py) to merge all epitopes.csv and all sources.csv into one file for every species.