# IEDB Protein Tree 

Mapping IEDB source antigens to genes and epitopes parent proteins. 


### Process
1. Fetch the species source antigens and epitopes from the IEDB MySQL backend.
2. Select the best proteome for that species from UniProt.
3. Assign genes to source antigens using BLAST and epitopes to their parent protein using PEPMatch.


### Inputs
- IEDB MySQL backend access
- List of IEDB species: [species.csv](species.csv)
    - This is updated with the [update_species.py](update_species.py) script
- `blastp` and `makeblastdb` binaries from [NCBI](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
- `hmmscan` built from [HMMER](http://hmmer.org/)
- Taxon ID for species to build for with `-t` flag or `-a` flag to build all species
- List of manual parents: [manual_assignments.csv](manual_assignments.csv)
    - This is the list of sources that have been manually assigned for their parents


### Running

To run the entire pipeline:

- for one species:
```bash
protein_tree/run_protein_tree.py -t <taxon ID>
```
- for all species:
```bash
protein_tree/run_protein_tree.py -a
```

each step can be run individually:

- fetch epitopes/source antigen data:
```bash
protein_tree/get_data.py -t <taxon ID>
```

- select best proteome:
```bash
protein_tree/select_proteome.py -t <taxon ID>
```

- assign gene and parent proteins:
```bash
protein_tree/assign_gene_protein.py -t <taxon ID>
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
- Add flag to do proteome selection - default is to use the proteome already selected (if there is one)
- Implement skipping a species if the data is the same as the last build
- Create manual assignments at the end of the all-species build
- Create a tree for visualization


### Workflow

<p align="center">
  <img src="docs/workflow.png">
</p>