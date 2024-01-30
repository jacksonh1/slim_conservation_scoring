- todo:
  - motif regex conservation scores
  - references
    - disorder matrix
    - iupred
    - capra/singh
  - generate "proteome" from orthodb in first step of pipeline
  - add optional link to odb database to get uniprot 2 odb id mapping
  - annotate msa's with organism names
  - add a parameter to specify which steps of the pipeline to run
  - Trim unused src code
  - slice by query gapless

# Table of Contents
- [Table of Contents](#table-of-contents)
- [motif conservation in disordered regions](#motif-conservation-in-disordered-regions)
- [setup TL;DR:](#setup-tldr)
- [database setup](#database-setup)
- [Using the pipeline](#using-the-pipeline)
  - [pipeline overview](#pipeline-overview)
  - [pipeline parameters](#pipeline-parameters)
    - [explanation of parameters](#explanation-of-parameters)
- [conservation scores](#conservation-scores)
  - [currently implemented scores:](#currently-implemented-scores)
    - [property entropy](#property-entropy)
    - [asymmetric sum-of-pairs score](#asymmetric-sum-of-pairs-score)
  - [multiple scores at once](#multiple-scores-at-once)
- [table annotations](#table-annotations)

# motif conservation in disordered regions

A series of tools to quantify the conservation of potential short linear motifs in disordered regions.<br>

- pipeline inputs:
  - A table containing candidate motifs and the name of the protein they are located in. 
    - The table can be generated in a number of ways, such as a regex search of a proteome or mapped sequences from a high-throughput experiment.
  - A "database" of multiple sequence alignments of each of the full length proteins in the input table and their orthologs.
- pipeline outputs:
  - a variety of conservation scores for each candidate motif, that are added back to the input table as a column

This repo includes an example database (`./data/example_orthogroup_database/`). The database is the output from my [orthoDB ortholog group preprocessing pipeline](https://github.com/jacksonh1/orthogroup_generation), and assumes that the ortholog groups were generated for all human proteins (in orthoDB) at different phylogenetic levels (Eukaryota, Metazoa, etc.) like in the [example script](https://github.com/jacksonh1/orthogroup_generation/tree/main/examples/ex3_all_human_genes). The basic idea is that you precalculate ortholog groups for conservation analysis and generate multiple sequence alignments for all of the groups. This pipeline is designed to access that "database" and calculate conservation scores for candidate motifs in a table. The pipeline can be used with any database of multiple sequence alignments, but the database key needs to be formatted in a specific way (see [database setup](#database-setup)).

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


## pipeline overview

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
The main pipeline is run via the file `./src/local_conservation_analysis_pipeline/conservation_pipeline.py`. The pipeline parameters are specified in a yaml file, for example:
```yaml
clear_files: true
output_folder: "./conservation_analysis"
database_filekey: "../../data/example_orthogroup_database/human_odb_groups/database_key.json"
table_file: "./table.csv"
hit_sequence_params:
  hit_sequence_search_method: "search"
idr_params:
  find_idrs: true
  # idr_map_file: "./idr_map.json"
  iupred_cutoff: 0.4
  gap_merge_threshold: 10
  idr_min_length: 8
filter_params:
  min_num_orthos: 20
new_score_methods:
  aln_asym_sum_of_pairs:
    matrix_name: "BLOSUM62_max_off_diagonal_norm"
    gap_frac_cutoff: 0.2
    overwrite: true
multilevel_plot_params:
  score_key: "aln_property_entropy"
  score_type: "zscore"
aln_slice_params:
  # n_flanking_aas: 20
  whole_idr: true
table_annotation_params:
  score_key_for_table: "property_entropy"
  motif_regex: "P.P.E"
  levels:
    - "Metazoa"
    - "Vertebrata"
  annotations:
    - "json_file"
    - "multi_level_plot"
    - "hit_start_position"
    - "regex"
    - "regex_match"
    - "regex_match_stpos_in_hit"
    - "conservation_string"
    - "aln_slice_file"
    # - "hit_scores"
    # - "hit_mean_score"
    # - "hit_z_scores"
    - "hit_mean_zscore"
    - "best mean z-score over 5 residue window"
    # - "best mean z-score over 10 residue window"
clean_analysis_files: false
```

### explanation of parameters
- `clear_files`: Whether or not to erase the `output_folder` (if it already exists) prior to running the pipeline
- `output_folder`: The folder where the analysis files will be saved. The pipeline creates a folder for each row in the input table, and saves the analysis files for that row in the folder.
- `database_filekey`: The file path to the json file that contains the database key. See [database setup](#database-setup) for more details.
- `table_file`: input table file, in csv format. The table needs to contain the following columns:
  - *gene_id* - the id of the full length protein. This is used to look up the database files in the `database_filekey`. It should correspond to the keys in the database key.
  - *hit_sequence* - the sequence of the candidate motif
  - *optional* - *hit start position*: the start position of the hit sequence in the full length protein. This is necessary when `hit_sequence_search_method` is set to `hit_sequence_positions`
  - *optional - hit end position*: the end position of the hit sequence in the full length protein
- `hit_sequence_params`:
  - `hit_sequence_search_method`: "search" or "given_positions" (default "search")
    - "search" - search for the hit sequence in the full length protein
    - "given_positions" - use the given start and end positions of the hit sequence in the full length protein
  - `longest_common_subsequence`: true or false (default true). If true, the longest common subsequence (lcs) of the hit sequence and the full length protein is used. If false, the hit sequence is used as is. Only used if `hit_sequence_search_method` is "search"
  - `lcs_min_length`: the minimum length of the longest common subsequence (default 20). If the lcs is shorter than this, the row is skipped. Only used if `longest_common_subsequence` is true and `hit_sequence_search_method` is "search"
  - `target_hit_length`: If the hit sequence is shorter than this, the pipeline will attempt to pad the hit sequence with flanking residues from the full length protein to achieve the `target_hit_length`. Only used if `hit_sequence_search_method` is "search" 
- `idr_params`:
  - `find_idrs`: true or false (default true). If true, the pipeline will use iupred to find the idrs in the full length protein. If false, the pipeline will use an external idr mapping file (`idr_map_file`).
  - `idr_map_file`: (not yet implemented). The file path to a json file that contains idr mappings. The file should have the gene_ids of the full length proteins as the keys, and the values should be a lists of lists, where each sublist contains the start and end positions of an idr. Required if `find_idrs` is false.
  - `iupred_cutoff`: the iupred cutoff to use when finding idrs (default 0.4). Scores above this are considered disordered.
  - `gap_merge_threshold`: if the distance between two idrs is less than or equal to this, the two idrs are merged into one (default 10)
  - `idr_min_length`: the minimum length of an idr (default 8). If an idr is shorter than this, it is discarded.
- `filter_params`:
  - `min_num_orthos`: minimum number of orthologs required for the group to be used in the analysis (default 20). For each phylogenetic level, if the corresponding alignment file has less than this number of orthologs, that level is skipped.
- `new_score_methods`:
  - Any new scores that are to be calculated are included here. The key is the score key corresponding to the conservation scoring method to use, and the value is a dictionary of parameters for the score. The parameters are specific to each score (See [scores](#conservation-scores)) but each score will have an input file, output file, and reference id which do **not** need to be provided here. See the [currently implemented scores](#currently-implemented-scores) for examples of how to format the yaml file.
- `multilevel_plot_params`: parameters for the multilevel plots of the hit sequence conservation
  - `score_key`: the score key to use for the multilevel plot (default "aln_property_entropy"). This can be a score that is already in the database key, or a new score that is calculated during pipeline execution.
  - `num_bg_scores_cutoff`: The minimum number or background scores required to calculate the zscores (default 20). If there are less than this number of background scores, the zscores are not calculated.
  - `score_type`: "score" or "zscore" (default "zscore"). If "score", the raw score is used. If "zscore", the zscore is used.
- `aln_slice_params`
  - `n_flanking_aas`: the number of residues to include on either side of the hit sequence in the output alignment slice file (default 20). It is the number of query sequence residues flanking the hit in the query sequence. So if there are gaps in the query sequence in the msa, those columns are not counted as flanking positions.
  - `whole_idr`: 
- `table_annotation_params`: parameters for adding conservation scores back to the input table
  - `score_key_for_table`: The score key corresponding to the score to add to the table (default "aln_property_entropy")
  - `motif_regex`: (not yet implemented). The regex to search for in the hit sequence (default None). If a regex is provided, an additional column is added to the table that that is the average conservation scoresof the residues in the hit sequence matching the regex. For example if the hit sequence is "PPPEQAPAPAEPGSA" and the regex is "P.P.E", the average conservation score of the residues "PAPAE" (xxxxxxPAPAExxxx) are calculated and added to the table.
  - `levels`: The phylogenetic levels to add to the table (default ["Metazoa", "Vertebrata"]). For each level, set of conservation scores is added to the table
  - `annotations`: The new columns to add to the output table. Some annotations that are available:
    - `json_file`: the path to the json file containing the conservation scores for the hit sequence
    - `multi_level_plot`: the path to the multilevel plot of the hit sequence conservation
    - `hit_start_position`: the start position of the hit sequence in the full length protein
    - `regex`: the regex used to find the motif within the the hit sequence
    - `regex_match`: the match of the regex in the hit sequence
    - `regex_match_stpos_in_hit`: the start position of the regex match within the hit sequence
    - `conservation_string`: a string of the conservation scores for the hit sequence. The string is the same length as the hit sequence. Residues with conservation z-scores < 0.5 are represented by a "_", and residues with z-scores >= 0.5 are represented by the corresponding amino acid.
    - `aln_slice_file`: the path to the alignment slice file
    - `hit_scores`: the conservation scores for the hit sequence (list that is the same length as hit sequence)
    - `hit_mean_score`: the mean conservation score for the hit sequence
    - `hit_z_scores`: the z-scores for the hit sequence (list that is the same length as hit sequence)
    - `hit_mean_zscore`: the mean z-score for the hit sequence
    - `best mean z-score over 5 residue window`: the best mean z-score over a 5 residue window
    - `best mean z-score over 10 residue window`: the best mean z-score over a 10 residue window
- `clean_analysis_files`: true or false (default false). If true, the analysis files are deleted after the pipeline is run (not yet implemented)

# conservation scores
The conservation score functions are located in `./src/local_conservation_scores/`. The main function for each score is in a separate file, and the functions are imported into a class in `./src/local_conservation_scores/__init__.py`. I am planning to add more scores in the future, and the plan is to add them to the class in `./src/local_conservation_scores/__init__.py` and create a new file for each score. I will include the docstring for each score below which will tell you what parameters are required. The docstrings are also available in the code of course.
<br>

## currently implemented scores:

### property entropy
- score from [capra and singh](https://pubmed.ncbi.nlm.nih.gov/17519246/)

docstring from `main` in `./src/local_conservation_scores/aln_property_entropy.py`:
```
calculate the property entropy score for each column in an alignment

Parameters
----------
input_alignment_file : str|Path
    input alignment file (fasta format)
output_file : str|Path
    output file (json). The output file will contain the following:
    "gap_mask" : list[bool]
        a list of bools indicating whether a column is masked by the gap mask
    "score_mask" : list[bool]
        a list of bools indicating whether a column is masked by the score mask. The score mask is the gap mask with additional columns masked if they are gaps in the reference sequence
    "scores" : list[float]
        a list of scores for each column in the alignment
reference_id : str
    the id of the sequence to use as the reference for the score mask (gaps in this sequence will be masked).
gap_frac_cutoff : float
    A number between 0 and 1. he fraction of gaps allowed in a column. If column has >gap_frac_cutoff gaps, it is masked in the gap mask
overwrite : bool, optional
    if True, overwrites the `output_file` if it already exists, by default False
```
example `new_score_methods` parameters in yaml file:
```yaml
new_score_methods:
  aln_property_entropy:
    gap_frac_cutoff: 0.2
    overwrite: true
```

### asymmetric sum-of-pairs score

- calculated using a substitution matrix. 
- It is asymmetric because only the "reference" sequence is compared to each ortholog (as opposed to an all-against-all comparison). The reference sequence is intended to be the protein (full length) that contains the hit sequence. The asymmetry makes it so that the conservation is specific to the reference sequence. For example, the `F` in this alignment is unique to the reference sequence, even though the position is generally conserved in the orthologs:
  ```
  reference   -------F--------
  ortholog 1  -------P--------
  ortholog 2  -------P--------
  ortholog 3  -------P--------
  ortholog 4  -------P--------
  ```
  The asymmetric score would give this position a low conservation score, whereas a symmetric score would give it a high conservation score. For this application we are interested in the conservation of the hit sequence in the reference sequence, so the asymmetric score is likely more appropriate.

docstring from `main` in `./src/local_conservation_scores/aln_asym_sum_of_pairs.py`:
```
calculate the asymmetric sum-of-pairs score for each column in an alignment

Parameters
----------
input_alignment_file : str|Path
    input alignment file (fasta format)
output_file : str|Path
    output file (json). The output file will contain the following:
    "gap_mask" : list[bool]
        a list of bools indicating whether a column is masked by the gap mask
    "score_mask" : list[bool]
        a list of bools indicating whether a column is masked by the score mask. The score mask is the gap mask with additional columns masked if they are gaps in the reference sequence
    "scores" : list[float]
        a list of scores for each column in the alignment
reference_id : str
    the id of the sequence to use as the reference for the score mask (gaps in this sequence will be masked).
matrix_name : str
    the name of the matrix to use for scoring.
    Available matrices:
        BLOSUM62_row_norm
        BLOSUM62_max_off_diagonal_norm
        EDSSMat50_row_norm
        EDSSMat50_max_off_diagonal_norm
gap_frac_cutoff : float
    A number between 0 and 1. he fraction of gaps allowed in a column. If column has >gap_frac_cutoff gaps, it is masked in the gap mask
overwrite : bool, optional
    if True, overwrites the `output_file` if it already exists, by default False
```

example `new_score_methods` parameters in yaml file:
```yaml
new_score_methods:
  aln_asym_sum_of_pairs:
    matrix_name: "BLOSUM62_max_off_diagonal_norm"
    gap_frac_cutoff: 0.2
    overwrite: true
```

## multiple scores at once
To calculate multiple conservation scores at once, you can add multiple score keys to the `new_score_methods` parameter. For example:
```yaml
new_score_methods:
  aln_asym_sum_of_pairs:
    matrix_name: "BLOSUM62_max_off_diagonal_norm"
    gap_frac_cutoff: 0.2
    overwrite: true
  aln_property_entropy:
    gap_frac_cutoff: 0.2
    overwrite: true
```




# table annotations
