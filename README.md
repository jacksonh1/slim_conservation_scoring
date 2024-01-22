# motif conservation in disordered regions

A series of tools to quantify the conservation of potential short linear motifs in disordered regions.<br>

- pipeline inputs:
  - A table containing candidate motifs and the name of the protein they are located in. 
    - The table can be generated in a number of ways, such as a regex search of a proteome or mapped sequences from a high-throughput experiment.
  - A "database" of multiple sequence alignments of each of the full length proteins in the input table and their orthologs.
- pipeline outputs:
  - a variety of conservation scores for each candidate motif, that are added back to the input table as a column

This repo includes an example database (`./data/example_orthogroup_database/`). The database is the output from my [orthoDB ortholog group preprocessing pipeline](https://github.com/jacksonh1/orthogroup_generation), and assumes that the ortholog groups were generated for all human proteins (in orthoDB) at different phylogenetic levels (Eukaryota, Metazoa, etc.) like in the [example script](https://github.com/jacksonh1/orthogroup_generation/tree/main/examples/ex3_all_human_genes).

# setup TL;DR:

1. download this repository
2. download iupred2a from [here](https://iupred2a.elte.hu/download_new)
3. edit the `.env` file with the location of the downloaded iupred2a folder: `ORTHODB_DATA_DIR=/absolute/path/to/iupred2a/folder/`
   - the folder should contain the file `iupred2a_lib.py`
4. create a new python environment with the dependencies: 
   - Mac - `conda env create -f environment.yml` <br>
   - Linux/windows WSL - `conda env create -f environment_linux.yml` <br>
5. activate the environment: `conda activate slim_conservation` <br>
6. install the src code as a local package: `pip install .` <br>

# database setup
An example database has already been prepared in this repo (`./data/example_orthogroup_database/human_odb_groups/`). The scripts used to generate this database are located in `./src/data_processing/`. The script `./src/data_processing/p1_initial_database_preparation/create_database.sh` was run with the [orthodb preprocessing pipeline](https://github.com/jacksonh1/orthogroup_generation), whereas the other scripts were run with this repo. They are numbered in the order they were run.

Regardless of whether or not you use the orthodb preprocessing pipeline, the database needs to consist of multiple sequence alignments for each of the full length proteins that are in the input table, and a json file (database key) that maps the gene id to the alignment files. The minimal json file should look like this:
```
{
    "gene_id": {
        "Vertebrata": {
            "alignment_file": "path/to/alignment/file/aln_file.fasta",
        },
        "Metazoa": {
            "alignment_file": "path/to/alignment/file/aln_file.fasta",
        },
        ...
    },
    ...
}
```
Where the gene_id corresponds to the id of the full length protein and the alignment file is a fasta file containing the multiple sequence alignment of the full length protein and its orthologs.

Important notes as of now:
- the pipeline assumes that the alignment is in fasta format.
- The gene_id in the database key is assumed to be the id of the full length protein in the fasta file.

Because the example database is created from the orthoDB, there are multiple phylogenetic levels for each gene (i.e. Vertebrata, Metazoa, etc.). If you are using a different source without levels, the json file should still be formatted the same way, but with just one "level". You can call the level whatever you like (such as "custom"). for example:
```
{
    "gene_id": {
        "custom": {
            "alignment_file": "path/to/alignment/file/aln_file.fasta",
        },
    },
    ...
}
```
For conservation scores that are calculated for an entire msa, the scores can be precalculated for all of the alignments in the database to save time when running the pipeline later. I did this for the example database with `./src/data_processing/p2_alignment_conservation_scores/property_entropy_scores.py`. Precalculated scores should be added to the database json file as a new key-value pair under the key `conservation_scores`. For example, I precalculated property entropy scores for the example database and the json file now looks like this:
```
{
    "gene_id": {
        "Vertebrata": {
            "alignment_file": "path/to/alignment/file/aln_file.fasta",
            "conservation_scores": {
                "property_entropy":"path/to/conservation/scores/file/conservation_scores.json",
            }
        },
        "Metazoa": {
            "alignment_file": "path/to/alignment/file/aln_file.fasta",
            "conservation_scores": {
                "property_entropy":"path/to/conservation/scores/file/conservation_scores.json",
            }
        },
        ...
    },
    ...
}
```
Whether or not to use a precalculated conservation score and any new scores to calculate are specified in the conservation pipeline parameters (see below).

# Using the pipeline
See `./examples/table_annotation/` for an example. <br>

The main pipeline is executed via the script `./`, which executes the following steps (code in `./src/local_conservation_analysis_pipeline/`) for each row in the input table:
1. setup the analysis folder (create folders and files for each row in the input table)
    - The gene id of the protein is looked up in the database key to find the alignment file(s) for the protein and its orthologs
2. define idr regions. 
    - By default this uses iupred (TODO: link) to define the idrs, however an external idr mapping file can be used instead.
3. find the hit sequence (candidate motif sequence) in the full length protein. This is done in one of two ways:
   - if you know the positions of the hit sequence in the full length protein, you can specify the start and end positions in the input table.
     - Be careful that the numbering actually lines up with the unaligned sequence in the alignment file. Different versions of the same protein can have slightly different numberings
   - search for an exact match in the full length sequence
     - if it is not found, the row is skipped
     - if it is found more than once, the row is skipped
     - In the future, I would like to add a parameter to add another row to the table if the hit sequence is found more than once
4. extract information from the alignment files
5. compute any new conservation scores that are to be calculated
6. construct conservation score plots of the hit sequence
7. create stylized html files of the alignment sliced around the hit sequence
8. add the conservation scores and file links back to the table and save a new annotated table.




## pipeline parameters
The pipeline parameters are specified in a yaml file




citations:
- disorder matrix
- iupred
- capra/singh
