# IEDB Protein Tree 

Mapping IEDB source antigens to genes and epitopes parent proteins. 

### Process
1. Fetch the species source antigens and epitopes from the IEDB MySQL backend.
2. Select the best proteome for that species from UniProt.
3. Assign genes to source antigens using BLAST and epitopes to their parent protein using PEPMatch.

### Inputs
- IEDB MySQL Backend access
- List of IEDB species: [species.csv](species.csv)
    - This is updated with the [update_species.py](update_species.py) script
- `blastp` and `makeblastdb` binaries from [NCBI](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
- Taxon ID for species to build for with `-t` flag or `-a` flag to build all species
- List of manual parents: [manual_assignments.csv](manual_assignments.csv)
    - This is the list of sources that have been manually assigned for their parents

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
- sources.csv - each source antigen with assigned gene
- epitopes.csv - each epitope with its source antigen and assigned parent protein
- [optional] gp_proteome.fasta - the gene priority proteome in FASTA if it exists

For all species:
- metrics.csv - the metadata from the build
    - Proteome ID
    - Proteome Taxon
    - Proteome Type
    - Source Antigen Count
    - Epitope Count
    - Successful Source Assignement (%)
    - Successful Epitope Assignment (%)
- all_epitopes.csv - combined epitope data
- all_sources.csv - combined source antigen data

Use [combine_data.py](combine_data.py) to merge all epitopes.csv and all sources.csv into one file for every species.

### TODO
- Improve BLAST assignments with protein existence levels / gene priority
- Isolate epitope matches to source protein assignment, not gene
- Create a tree for visualization
- Add higher taxon ranks for the final tree