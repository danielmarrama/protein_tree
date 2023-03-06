# IEDB Protein Tree 

Mapping epitopes and source antigens to their parent protein (and gene).

### Process
1. Fetch epitopes and source antigens for each species in the IEDB.
2. Select the best proteome for that species.
3. Assign genes to source antigens using BLAST and PEPMatch.
4. Assign epitopes to their parent protein. 

### Inputs
- List of IEDB species: [species.csv](species.csv)
- Taxon ID for species to build or `all` to build all species. 

### Running
For one species:
``` python
./protein_tree.py -u <IEDB backend username> -p <IEDB backend password> -t <taxon ID>
```
or for all species:
``` python
./protein_tree.py -u <IEDB backend username> -p <IEDB backend password> -a
```

### Outputs

Mappings:
Epitope --> Source Antigen
Source Antigen --> Gene
Gene --> Isoforms (Proteins)
Epitope --> Parent Protein