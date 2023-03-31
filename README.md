# IEDB Protein Tree 

Mapping source antigens to genes and epitopes to parent proteins.

### Process
1. Fetch epitopes and source antigens for species in the IEDB.
2. Select the best proteome for that species.
3. Assign genes to source antigens using BLAST and epitopes to their parent protein using PEPMatch.

### Inputs
- List of IEDB species: [species.csv](species.csv)
- IEDB MySQL Backend username and password
- `blastp` and `makeblastdb` executables from [NCBI](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
- Taxon ID for species to build for with `-t` flag or `-a` flag to build all species.
- List of manual parents: [manual_assignments.csv](manual_assignments.csv)
    - This is the list of sources that have been manually assigned for their parents.

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
- [optional] no_blast_match_ids.txt - list of source IDs that did not get any BLAST matches
- [optional] sources_missing_seqs.txt - sources from source table that did not have a sequence

For all species:
- metrics.csv - the metadata from the build
    - Proteome ID
    - Proteome Taxon
    - Proteome Type
    - \# of source antigens
    - \# of epitopes
    - \# of sources missing sequences
    - \# of sources with no BLAST match
    - \# of sources with a BLAST match
    - \# of epitopes with a PEPMatch match
    - % of successful gene assignments
    - % of successful epitope assignments
- all_epitopes.csv - combined epitope data
- all_sources.csv - combined source antigen data

Use [combine_data.py](combine_data.py) to merge all epitopes.csv and all sources.csv into one file for every species.

### TODO
- Make a mapping of all children to parent species
- Fix epitope data pull
- Use manual_assignments.csv to override parent and gene assignments
- Add step to handle allergens using labels from IUIS
- Add step to handle MHC molecules from MRO
- Use PEPMatch to search discontinuous epitopes
- Test other alignment tools (MMseqs2 or DIAMOND) for accuracy
