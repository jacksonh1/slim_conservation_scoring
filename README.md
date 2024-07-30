Written by Jackson Halpin <br>

This work was supported by the National Institutes of Health under Award Number R35GM149227. The content is solely the responsibility of the authors and does not necessarily represent the official views of the National Institutes of Health.

---

# Table of Contents
- [Table of Contents](#table-of-contents)
- [motif conservation in disordered regions](#motif-conservation-in-disordered-regions)
- [note on pipeline structure](#note-on-pipeline-structure)
- [setup TL;DR:](#setup-tldr)
- [database setup](#database-setup)
- [Using the pipeline](#using-the-pipeline)
  - [pipeline overview](#pipeline-overview)
  - [pipeline parameters](#pipeline-parameters)
    - [explanation of parameters](#explanation-of-parameters)
- [conservation scores](#conservation-scores)
  - [currently implemented MSA scores:](#currently-implemented-msa-scores)
    - [shannon entropy](#shannon-entropy)
    - [property entropy](#property-entropy)
    - [asymmetric sum-of-pairs score](#asymmetric-sum-of-pairs-score)
  - [multiple scores](#multiple-scores)
  - [currently implemented pairk alignment methods:](#currently-implemented-pairk-alignment-methods)
    - [pairk alignment](#pairk-alignment)
    - [pairk alignment with needleman-wunsch](#pairk-alignment-with-needleman-wunsch)
    - [pairk alignment with embeddings](#pairk-alignment-with-embeddings)
  - [multiple pairk alignment methods](#multiple-pairk-alignment-methods)
  - [currently implemented pairk conservation methods:](#currently-implemented-pairk-conservation-methods)
    - [pairk conservation](#pairk-conservation)
  - [references](#references)


# motif conservation in disordered regions

A series of tools to quantify the conservation of potential short linear motifs in disordered regions.<br>

The basic idea is that you precalculate ortholog groups for conservation analysis and generate multiple sequence alignments for all of the groups. This pipeline is designed to access that "database" and calculate conservation scores for candidate motifs in a table. The pipeline can be used with any database of multiple sequence alignments, but the database key needs to be formatted in a specific way (see [database setup](#database-setup)). This repo includes an example database (`./data/example_orthogroup_database/`). The database is the output from my [orthoDB ortholog group preprocessing pipeline](https://github.com/jacksonh1/slim_conservation_orthogroup_generation), and assumes that the ortholog groups were generated for all human proteins (in orthoDB) at different phylogenetic levels (Eukaryota, Metazoa, etc.) like in the [example](https://github.com/jacksonh1/slim_conservation_orthogroup_generation/tree/main/examples/ex3_all_human_genes). The code for generating the example database is in `./data_processing/`.


- pipeline inputs:
  - A table containing candidate motifs and the name of the protein they are located in. 
    - The table can be generated in a number of ways, such as a regex search of a proteome or mapped sequences from a high-throughput experiment.
  - A "database" of multiple sequence alignments of each of the full length proteins in the input table and their orthologs.
- pipeline outputs:
  - a variety of conservation scores for each candidate motif, that are added back to the input table as columns

conservation can be calculated via traditional MSA-based methods or via [pairk](https://github.com/jacksonh1/pairk).
See [conservation scores](#conservation-scores) for more information on the conservation scores that can be calculated and see [pairk's documentation](https://pairk.readthedocs.io/en/latest/) for more information on pairk.

# note on pipeline structure
- I wanted to avoid using a workflow manager (e.g. snakemake or nextflow) and a more sophisticated database structure because I wanted to keep the pipeline easy for anyone to use, however, I think the result is potentially more confusing. If I had the time to redo it, I would probably use a workflow manager to handle the pipeline and a different database structure to store the data and results.

# setup TL;DR:

1. download this repository
2. download iupred2a from [here](https://iupred2a.elte.hu/download_new)
3. edit the `.env` (`./slim_conservation_scoring/env_variables/.env`) file with the location of the downloaded iupred2a folder: `IUPRED2A_LIB_DIR=/absolute/path/to/iupred2a/folder/`
   - the folder should contain the file `iupred2a_lib.py`
4. create a new python environment with the dependencies: 
   - Mac - `conda env create -f environment_mac.yml` <br>
   - Linux/windows WSL - `conda env create -f environment_linux.yml` <br>
5. activate the environment: `conda activate slim_conservation_scoring` <br>
6. install the src code as a local package: `pip install .` <br>

# database setup
An example database has already been prepared in this repo (`./data/example_orthogroup_database/human_odb_groups/`). The scripts used to generate this database are located in `./data_processing/` and were run with the [orthodb orthogroup preprocessing pipeline](https://github.com/jacksonh1/orthogroup_generation).

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
<!-- For conservation scores that are calculated for an entire msa, the scores can be precalculated for all of the alignments in the database to save time when running the pipeline later. I did this for the example database with `./src/data_processing/p2_alignment_conservation_scores/property_entropy_scores.py`. Precalculated scores should be added to the database json file as new entries under the key `conservation_scores`. Each score requires a key and associated "file", "score_function_name", and "score_params" ent. For example, I precalculated property entropy scores for the example database and the json file now looks like this:
```
{
    "gene_id": {
        "Vertebrata": {
            "alignment_file": "path/to/alignment/file/aln_file.fasta",
            "conservation_scores": {
                "aln_property_entropy": {
                    "file": "path/to/conservation/scores/file/conservation_scores.json",
                    "score_function_name": "aln_property_entropy",
                    "score_params": {
                        "gap_frac_cutoff": 0.2,
                        "overwrite": true
                    }
                }
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
The precalculated conservation score can just be specified in the conservation pipeline parameters (see below). -->

# Using the pipeline
See `./examples/table_annotation/` for an example. <br>


## pipeline overview

The main pipeline is executed via the script `./slim_conservation_scoring/pipeline/conservation_pipeline.py`, which executes the following steps (code in `./slim_conservation_scoring/pipeline/`) for each row in the input table:
- s1. setup the analysis folder (create folders and files for each row in the input table)
    - a reference index is used to keep track of each row in the table. Each row is associated with a unique reference index.
    - For each row, the gene id of the protein is looked up in the database key to find the alignment file(s) for the protein and its orthologs
    - a folder is created for the row in the output folder, and a json file is created in the folder that stores all of the information for the row. The json file is used to keep track of the analysis files for the row.
- s2. define idr regions. 
    - By default this uses iupred (TODO: link) to define the idrs, however an external idr mapping file can be used instead.
    - the json file is updated with the IDR information
- s3. find the hit sequence (candidate motif sequence) in the full length protein. This is done in one of two ways:
   - if you know the positions of the hit sequence in the full length protein, you can specify the start and end positions in the input table.
     - Be careful that the numbering actually lines up with the unaligned sequence in the alignment file. Different versions of the same protein can have slightly different numberings
   - search for an exact match in the full length sequence
     - if it is not found, the row is skipped
     - if it is found more than once, the row is skipped
     - (In the future, I would like to add a parameter to add another row to the table if the hit sequence is found more than once)
   - the json file is updated with the hit sequence information and whether or not the hit is contained within an IDR
- s4. extract information from the alignment files and database and add it to the json file for each row.
- s5. compute any new conservation scores that are to be calculated. Add the scores to the json file for each row.
- s6. construct conservation score plots of the hit sequence
- s7. create stylized html files of the alignment sliced around the hit sequence
- s8. calculate table annotations to be added back to the input table
- s9. add the conservation scores and file paths back to the table and save a new annotated table.


## pipeline parameters
The main pipeline is run via the file `./slim_conservation_scoring/pipeline/conservation_pipeline.py`. The pipeline parameters are specified in a yaml file, for example (example is default configuration; all entries except database_filekey, table_file, and output_folder are optional):
```yaml
database_filekey: "../../data/example_orthogroup_database/human_odb_groups/database_key.json"
table_file: "./table.csv"
output_folder: "./conservation_analysis"
clear_files: True
new_index: True
steps_to_run:
  - "s1" # setup folders
  - "s2" # define_idrs
  - "s3" # find_hit
  - "s4" # add_lvlinfo
  - "s5" # scores
  - "s6" # multilevel_plots
  - "s7" # output_aln_slice
  - "s8" # calculate_annotations
  - "s9" # add_annotations2table
hit_sequence_params:
  hit_sequence_search_method: "search"
  longest_common_subsequence: false
  lcs_min_length: 20
  target_hit_length: 0
idr_params:
  find_idrs: true
  idr_map_file: null
  iupred_cutoff: 0.4
  gap_merge_threshold: 10
  idr_min_length: 8
filter_params:
  min_num_orthos: 20
msa_score_methods:
  aln_property_entropy:
    function_name: "aln_property_entropy"
    function_params:
      gap_frac_cutoff: 0.2
      overwrite: false
  aln_shannon_entropy:
    function_name: "aln_shannon_entropy"
    function_params:
      gap_frac_cutoff: 0.2
      overwrite: false
pairk_aln_methods:
  pairk_aln_lf5_rf5_edssmat50:
    function_name: "pairk_aln"
    function_params:
      overwrite: false
      matrix_name: "EDSSMat50"
    lflank: 5
    rflank: 5
esm_params:
  processes: 4
  threads: 1
  device: "cuda"
  model_name: "esm2_t33_650M_UR50D"
pairk_embedding_aln_methods:
  pairk_aln_embedding_lf5_rf5:
    function_name: "pairk_aln_embedding"
    function_params:
      overwrite: false
    lflank: 5
    rflank: 5
    level: "Vertebrata"
pairk_conservation_params:
  kmer_conservation_function_name: "pairk_conservation"
  columnwise_score_function_name: "shannon_entropy"
  bg_cutoff: 50
  bg_kmer_cutoff: 10
multilevel_plot_params:
  score_key: "aln_property_entropy"
  num_bg_scores_cutoff: 20
  score_type: "z_score"
  strip_gaps: false
aln_slice_params:
  n_flanking_aas: 20
  whole_idr: false
table_annotation_params:
  score_key_for_table: "aln_property_entropy"
  motif_regex: "P.P.E"
  levels:
    - "Metazoa"
    - "Vertebrata"
  annotations: 
    - 'json_file'
clean_analysis_files: false
```

### explanation of parameters
- `steps_to_run`: The steps of the pipeline to execute (corresponding to the steps in [pipeline overview](#pipeline-overview)). It will run all steps by default.
- `clear_files`: Whether or not to erase the `output_folder` (if it already exists) prior to running the pipeline
- `new_index`: Whether or not to add a new reference index column to the table. The reference index is used by the pipeline to keep track of the rows in the table. 
  - If `new_index`=true and `clear_files`=true, the table is copied into the `output_folder` (with the suffix `_original_reindexed` added to the filename) and a new reference index column is added to the copied table.
  - If `new_index`=true and `clear_files`=false, an error is raised to avoid really confusing errors. A new reference index column cannot be created without making sure there are no files in `output_folder` first.
  - If `new_index`=false, the pipeline will first look for a reindexed table file in the `output_folder`. 
    - If the reindexed table file exists, it will use the reindexed table file and it's reference index column.
    - If the reindexed table does not exist, the pipeline will assume that the input table already has a reference_index column (must be named **reference_index**) and use that column as the reference index without reindexing it. It will still created a copy of this table and save it in the `output_folder` with the suffix `_original_reindexed` added to the filename.
- `output_folder`: The folder where the analysis files will be saved. If `clear_files`=true, this entire folder is deleted. The pipeline creates an additional folder within this folder for each row in the input table, and saves the analysis files for that row in the folder. 
- `database_filekey`: The file path to the json file that contains the database key. See [database setup](#database-setup) for more details.
- `table_file`: input table file, in csv format. The table needs to contain the following columns:
  - *gene_id* - the id of the full length protein. This is used to look up the database files in the `database_filekey`. It should correspond to the keys in the database key.
  - *hit_sequence* - the sequence of the candidate motif
  - (optional) *hit start position*: the start position of the hit sequence in the full length protein. This is necessary when `hit_sequence_search_method` is set to "given_positions"
  - (optional) *hit end position*: the end position of the hit sequence in the full length protein. This is necessary when `hit_sequence_search_method` is set to "given_positions"
- `idr_params`:
  - `find_idrs`: true or false (default true). If true, the pipeline will use iupred to find the idrs in the full length protein. If false, the pipeline will use an external idr mapping file (`idr_map_file`).
  - `idr_map_file`: (not yet implemented). The file path to a json file that contains idr mappings. The file should have the gene_ids of the full length proteins as the keys, and the values should be a lists of lists, where each sublist contains the start and end positions of an idr. Required if `find_idrs` is false.
  - `iupred_cutoff`: the iupred cutoff to use when finding idrs (default 0.4). Scores above this are considered disordered.
  - `gap_merge_threshold`: if the distance between two idrs is less than or equal to this, the two idrs are merged into one (default 10)
  - `idr_min_length`: the minimum length of an idr (default 8). If an idr is shorter than this, it is discarded.
- `filter_params`:
  - `min_num_orthos`: minimum number of orthologs required for the group to be used in the analysis (default 20). For each phylogenetic level, if the corresponding alignment file has less than this number of orthologs, that level is skipped.
- `hit_sequence_params`:
  - `hit_sequence_search_method`: "search" or "given_positions" (default "search")
    - "search" - search for the hit sequence in the full length protein
    - "given_positions" - uses given start and end positions of the hit sequence (in the full length protein) provided by the user as columns in the input table.
  - `longest_common_subsequence`: true or false (default true). If true, the longest common subsequence (lcs) of the hit sequence and the full length protein is used. If false, the hit sequence is used as is. Only used if `hit_sequence_search_method` is "search"
  - `lcs_min_length`: the minimum length of the longest common subsequence (default 20). If the lcs is shorter than this, the row is skipped. Only used if `longest_common_subsequence` is true and `hit_sequence_search_method` is "search"
  - `target_hit_length`: If the hit sequence is shorter than this, the pipeline will attempt to pad the hit sequence with flanking residues from the full length protein to achieve the `target_hit_length`. Only used if `hit_sequence_search_method` is "search" 
- `msa_score_methods`:
  - Any new msa-based scores that are to be calculated are included here. For each scoring method to be calculated, the scoring method is provided in the following format:
      ```yaml
      msa_score_methods:
        score_key1:
          function_name: "function_name"
          function_params:
            param1: value1
            param2: value2
        score_key2:
          function_name: "function_name"
          function_params:
            param1: value1
            param2: value2
          level: "Vertebrata"
      ```
  - `score_key` (score_key1/score_key2): A unique identifier for each individual score (including pairk scores). can be any string but must be unique.
    - `function_name`: The name of the conservation scoring function to use (See [scores](#conservation-scores))
    - `function_params`: Keyword arguments specific to each score function (See [scores](#conservation-scores)). Each scoring method will have an input file, output file, and reference id which do **not** need to be provided here. 
    - `level`: The level in the database to use for the score calculation. If not provided, the pipeline will calculate the score for every level in the database.
  - See the [currently implemented scores](#currently-implemented-scores) for examples of how to format the yaml file.
- `pairk_aln_methods`:
  - Here, the parameters for [pairk](https://github.com/jacksonh1/pairk)'s k-mer alignment step are specified. The parameters are entered in a similar format as the `msa_score_methods`:
      ```yaml
      pairk_aln_methods:
        score_key3:
          function_name: "function_name"
          function_params:
            param1: value1
            param2: value2
          lflank: 5
          rflank: 5
        score_key4:
          function_name: "function_name"
          function_params:
            param1: value1
            param2: value2
          lflank: 0
          rflank: 0
          level: "Vertebrata"
      ```
  - `score_key` (score_key3/score_key4): A unique identifier for each individual score. can be any string but must be unique.
    - `function_name`: The name of the conservation scoring function to use (See [scores](#conservation-scores))
    - `function_params`: Keyword arguments specific to each score function (See [scores](#conservation-scores)). Each scoring method will have an input file, output file, and reference id which do **not** need to be provided here.
    - `lflank` and `rflank`: the number of additional residues (from the full length sequence) to include on the left (`lflank`) and right (`rflank`) of the hit sequence. This will change the value of k in the pairk k-mer analysis. 
    - `level`: The level in the database to use for the score calculation. If not provided, the pipeline will calculate the score for every level in the database.
- `pairk_embedding_aln_methods`:
  - Any embedding-based `pairk` alignment methods are added here. The parameters are formated in the same way as the `pairk_aln_methods`
- `esm_params`:
  - `processes`: the number of processes to use for the ESM embedding (default 4). This is the number of concurrent table rows that are processed at once with the embedding method. It is included here as a separate parameter because the ESM embedding method uses a lot of memory. This parameter overrides the `n_cores` parameter in the main pipeline function (`./slim_conservation_scoring/pipeline/conservation_pipeline.py`) and command line script (`./slim_conservation_scoring/scripts/conservation_analysis.py`).
  - `threads`: the number of threads to use for the ESM embedding (default 1)
  - `device`: the device to use for the ESM embedding, must be 'cpu' or 'cuda' (default "cuda")
- `pairk_conservation_params`: parameters for calculating conservation scores from the pairk alignment results. Only applies to results from `pairk_aln_methods` and `pairk_embedding_aln_methods`.
  - `kmer_conservation_function_name`: the name of the function to use for calculating the conservation of the k-mers in the pairk alignment (default "pairk_conservation"). see [scores](#conservation-scores) for more information.
  - `columnwise_score_function_name`: the name of the function to use for calculating the conservation of the columns in the pairk alignment (default "shannon_entropy"). see [scores](#conservation-scores) for more information.
  - `bg_cutoff`: the minimum number of background scores required to calculate the z-scores (default 50). If there are less than this number of background scores, the z-scores are not calculated.
  - `bg_kmer_cutoff`: the minimum number of background kmers required to calculate the z-scores (default 10). If there are less than this number of background kmers, the z-scores are not calculated.
- `multilevel_plot_params`: parameters for the multilevel plots of the hit sequence conservation.
  - `score_key`: the score key to use for the multilevel plot (default "aln_property_entropy"). Must be a score_key that is specified via `msa_score_methods`, `pairk_aln_methods`, or `pairk_embedding_aln_methods`.
  - `num_bg_scores_cutoff`: The minimum number or background scores required to calculate the z-scores (default 20). If there are less than this number of background scores, the z-scores are not calculated.
  - `score_type`: "score" or "z_score" (default "z_score"). If "score", the raw score is used. If "z_score", the z_score is used.
- `aln_slice_params`:
  - `n_flanking_aas`: the number of residues to include on either side of the hit sequence in the output alignment slice file (default 20). It is the number of query sequence residues flanking the hit in the query sequence. So if there are gaps in the query sequence in the msa, those columns are not counted as flanking positions.
  - `whole_idr`: whether or not to include the entire IDR in the alignment slice: Either true or false (default false).
- `table_annotation_params`: parameters for adding conservation scores back to the input table.
  - `score_key_for_table`: The score key corresponding to the score to add to the table (default "aln_property_entropy"). Must be a score_key that is specified via `msa_score_methods`, `pairk_aln_methods`, or `pairk_embedding_aln_methods`.
  - `motif_regex`: The regex to search for in the hit sequence (default None). If a regex is provided, additional columns can be added to the table via the regex options below. 
  - `levels`: The phylogenetic levels to add to the table (default ["Metazoa", "Vertebrata"]). For each level, the pipeline will add annotations to the table.
  - `annotations`: The new columns to add to the output table. Some annotations are calculated for each "level" specified above. Available annotations:
    - `json_file`: the path to the json file containing the conservation scores for the hit sequence
    - `multi_level_plot`: the path to the multilevel plot of the hit sequence conservation
    - `hit_start_position`: the start position of the hit sequence in the full length protein
    - `regex`: the regex used to find the motif within the the hit sequence
    - `regex_match`: the match of the regex in the hit sequence. i.e. the sequence within the hit sequence matching the regex
    - `regex_match_stpos_in_hit`: the start position of the regex match within the hit sequence
    - `regex_match_scores`: the conservation scores for the residues in the regex match (list)
    - `regex_match_mean_score`: the mean conservation score of the residues in the regex match
    - `regex_match_z_scores`: the z-scores for the residues in the regex match (list)
    - `regex_match_mean_zscore`: the mean z-score of the residues in the regex match
    - `conservation_string`: a string of the conservation scores for the hit sequence. The string is the same length as the hit sequence. Residues with conservation z-scores < 0.5 are represented by a "_", and residues with z-scores >= 0.5 are represented by the corresponding amino acid.
    - `aln_slice_file`: the path to the alignment slice file
    - `hit_scores`: the conservation scores for the hit sequence (list that is the same length as hit sequence)
    - `hit_mean_score`: the mean conservation score for the hit sequence
    - `hit_z_scores`: the z-scores for the hit sequence (list that is the same length as hit sequence)
    - `hit_mean_zscore`: the mean z-score for the hit sequence
    - `best mean z-score over 5 residue window`: the best mean z-score over a 5 residue window
    - `best mean z-score over 10 residue window`: the best mean z-score over a 10 residue window
- `clean_analysis_files`: true or false (default false). If true, the analysis files are deleted after the pipeline is run


# conservation scores
The conservation score functions are located in `./slim_conservation_scoring/conservation_scores/`. The main function for each score is in a separate file, and the functions are imported into classes in `./slim_conservation_scoring/conservation_scores/__init__.py`. The `function_name` parameter under `msa_score_methods`, `pairk_aln_methods`, or `pairk_embedding_aln_methods` specifies the function to run (it is used to access the function from the class in `./slim_conservation_scoring/conservation_scores/__init__.py`). Similarly, the `kmer_conservation_function_name` and the `columnwise_score_function_name` parameters under `pairk_conservation_params` also specify functions to use for the pairk conservation calculations.
<br>

classes in `./slim_conservation_scoring/conservation_scores/__init__.py` and their corresponding parameters in the configuration file:
- `MSAScoreMethods`: class that contains the functions for calculating conservation scores from multiple sequence alignments. 
  - available functions:
    - `aln_asym_sum_of_pairs`
    - `aln_property_entropy`
    - `aln_shannon_entropy`
  - **Yaml file** - `function_name` parameter under `msa_score_methods`. `function_params` are passed to the function as keyword arguments.
- `PairKmerAlnMethods`: class that contains the functions for running the pairk k-mer alignment. 
  - available functions:
    - `pairk_aln`
    - `pairk_aln_needleman`
    - `pairk_aln_embedding`
  - **Yaml file** - `function_name` parameter under `pairk_aln_methods` and `pairk_embedding_aln_methods`. `function_params` are passed to the function as keyword arguments.
- `PairKmerConservationMethods`: class that contains the functions for calculating conservation scores from the pairk k-mer alignment. There is currently only 1 function in this class.
  - available functions:
    - `pairk_conservation`
  - **Yaml file** - `kmer_conservation_function_name` parameter under `pairk_conservation_params`. `bg_cutoff` and `bg_kmer_cutoff` are used as keyword arguments.
- `ColumnwiseScoreMethods`: class that contains functions that calculate a conservation score from a sequence alignment column (used as a parameter for pairk's k-mer conservation methods). 
  - available functions:
    - `shannon_entropy`
    - `property_entropy`
  - **Yaml file** - `columnwise_score_function_name` parameter under `pairk_conservation_params`.


<br>

## currently implemented MSA scores:
For conservation scores calculated from MSAs, the background scores for the z-score calculation are the conservation scores of every column in the IDR region of the MSA that is not masked. An MSA column is masked if it corresponds to a gap in the query sequence or if the column has more than the allowed fraction of gaps (>`gap_frac_cutoff`, parameter in the below functions). 

### shannon entropy
- score (including code) is from [capra and singh](https://pubmed.ncbi.nlm.nih.gov/17519246/)

docstring from `main` in `./slim_conservation_scoring/conservation_scores/aln_shannon_entropy.py`:
```
calculate the shannon entropy score for each column in an alignment

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
**`reference_id`, `output_file`, and `input_alignment_file` should not be specified in the yaml file. They are handled by the pipeline. other keyword arguments can be specified in the yaml file under the `function_params` keyword**

example `msa_score_methods` parameters in yaml file:
```yaml
msa_score_methods:
  msa_SE:
    function_name: "aln_shannon_entropy"
    function_params: 
      gap_frac_cutoff: 0.2
      overwrite: false
```

### property entropy
- score (including code) is from [capra and singh](https://pubmed.ncbi.nlm.nih.gov/17519246/)

docstring from `main` in `./slim_conservation_scoring/conservation_scores/aln_property_entropy.py`:
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
**`reference_id`, `output_file`, and `input_alignment_file` should not be specified in the yaml file. They are handled by the pipeline. other keyword arguments can be specified in the yaml file under the `function_params` keyword**

example `msa_score_methods` parameters in yaml file:
```yaml
msa_score_methods:
  msa_PE:
    function_name: "aln_property_entropy"
    function_params: 
      gap_frac_cutoff: 0.2
      overwrite: false
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

docstring from `main` in `./slim_conservation_scoring/conservation_scores/aln_asym_sum_of_pairs.py`:
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
**`reference_id`, `output_file`, and `input_alignment_file` should not be specified in the yaml file. They are handled by the pipeline. other keyword arguments can be specified in the yaml file under the `function_params` keyword**

example `msa_score_methods` parameters in yaml file:
```yaml
msa_score_methods:
  msa_sop:
    function_name: "aln_asym_sum_of_pairs"
    function_params: 
      gap_frac_cutoff: 0.2
      overwrite: false
      matrix_name: "BLOSUM62_max_off_diagonal_norm"
```
## multiple scores
To calculate multiple conservation scores at once, you can add multiple score keys to the `msa_score_methods` parameter. For example:
```yaml
msa_score_methods:
  msa_sop:
    function_name: "aln_asym_sum_of_pairs"
    function_params: 
      gap_frac_cutoff: 0.2
      overwrite: false
      matrix_name: "BLOSUM62_max_off_diagonal_norm"
  msa_PE:
    function_name: "aln_property_entropy"
    function_params: 
      gap_frac_cutoff: 0.2
      overwrite: false
```

## currently implemented pairk alignment methods:
There are different methods to run the k-mer alignment step of pairk. See [pairk](https://github.com/jacksonh1/pairk) for more details about pairk alignment.

### pairk alignment

docstring from `main` in `./slim_conservation_scoring/conservation_scores/pairk_aln.py`:
```
run pairk alignment from a fasta alignment file

Parameters
----------
input_alignment_file : str | Path
    input fasta alignment file
output_file : str | Path
    output file for pairk alignment results in json format
reference_id : str
    the sequence id of the query sequence
k : int
    the value of k for pairk alignment
idr_aln_st : int
    start position of the IDR in the alignment
idr_aln_end : int
    end position of the IDR in the alignment
overwrite : bool, optional
    whether or not to overwrite the output file if it exists, by default False
matrix_name : str, optional
    the name of the scoring matrix to use for pairk.pairk_alignment, by default "EDSSMat50". The available matrices can be viewed by importing pairk and running the function `pairk.print_available_matrices()`.
```
**`k`, `reference_id`, `output_file`, and `input_alignment_file` should not be specified in the yaml file. They are handled by the pipeline. other keyword arguments can be specified in the yaml file under the `function_params` keyword**

example `pairk_aln_methods` parameters in yaml file:
```yaml
pairk_aln_methods:
  pairk_aln_lf5_rf5_blos62:
    function_name: "pairk_aln"
    function_params: 
      overwrite: false
      matrix_name: "BLOSUM62"
    lflank: 5
    rflank: 5
    level: "Vertebrata"
```

### pairk alignment with needleman-wunsch
This is the same as the pairk alignment, but uses the needleman-wunsch algorithm instead of the default exhaustive algorithm. The yaml file parameters are the same as the pairk alignment, except that the `function_name` is "pairk_aln_needleman".

### pairk alignment with embeddings
This function has the same parameters as [pairk alignment](#pairk-alignment), except that it does not take a `matrix_name` parameter. Sequence embeddings are calculated using the ESM model, and the embeddings are used to run the pairk alignment. The embeddings are calculated for the query sequence and the orthologs, and the pairk alignment is calculated using the embeddings.

example `pairk_embedding_aln_methods` parameters in yaml file:
```yaml
pairk_embedding_aln_methods:
  pairk_aln_embedding_lf0_rf0:
    function_name: "pairk_aln_embedding"
    function_params:
      overwrite: false
    lflank: 0
    rflank: 0
    level: "Vertebrata"
```

## multiple pairk alignment methods
To run multiple pairk alignment methods at once, you can add multiple score keys to the `pairk_aln_methods` or `pairk_embedding_aln_methods` parameter. For example:
```yaml
pairk_aln_methods:
  pairk_aln_lf5_rf5_edssmat50:
    function_name: "pairk_aln"
    function_params:
      overwrite: false
      matrix_name: "EDSSMat50"
    lflank: 5
    rflank: 5
    level: "Vertebrata"
  pairk_aln_lf0_rf0_edssmat50:
    function_name: "pairk_aln"
    function_params:
      overwrite: false
      matrix_name: "EDSSMat50"
    lflank: 0
    rflank: 0
pairk_embedding_aln_methods:
  pairk_aln_embedding_lf5_rf5:
    function_name: "pairk_aln_embedding"
    function_params:
      overwrite: false
    lflank: 5
    rflank: 5
  pairk_aln_embedding_lf0_rf0:
    function_name: "pairk_aln_embedding"
    function_params:
      overwrite: false
    lflank: 0
    rflank: 0
    level: "Vertebrata"
```

## currently implemented pairk conservation methods:
k-mer conservation is calculated from the results of the pairk alignment step. There are some additional parameters that can be specified here, so I did not combine pairk alignment and pairk conservation into 1 step. The pairk conservation is configured via the `pairk_conservation_params` key of the config file. See [pairk](https://github.com/jacksonh1/pairk) for more details about pairk conservation.

### pairk conservation
There is only one currently implemented function and it's likely to stay that way.

docstring from `pairk_conservation_from_json` in `./slim_conservation_scoring/conservation_scores/pairk_conservation.py`:
```
Calculate conservation scores from the pairk alignment results using the `pairk.calculate_conservation` function. Only the scores from the hit k-mer are returned.

Parameters
----------
kmer_aln_json : str | Path
    The json file storing the pairk alignment results
hit_position : int
    The position of the hit k-mer in the query IDR sequence
columnwise_score_func : Callable, optional
    A function to calculate conservation scores in a columnwise manner, can
    be any function that takes a string representing a column of an MSA.
    by default it is the property_entropy function from Capra and Singh 2007,
    DOI: 10.1093/bioinformatics/btm270.
bg_cutoff : int, optional
    the minimum number of background scores required to calculate the
    z-scores, by default 50
bg_kmer_cutoff : int, optional
    the minimum number of background kmers required to calculate the
    z-scores, by default 10

Returns
-------
PairwiseScoreResults
    Object containing pairk conservation results for the hit k-mer.
    Includes the hit k-mer sequence, the conservation scores, the conservation
    z-scores, and the background scores.

Raises
------
ValueError
    If there are fewer than `bg_kmer_cutoff` kmers to use for background scores
ValueError
    If there are fewer than `bg_cutoff` background scores
```
example `pairk_conservation_params` parameters in yaml file:
```yaml
pairk_conservation_params:
  kmer_conservation_function_name: "pairk_conservation"
  columnwise_score_function_name: "shannon_entropy"
  bg_cutoff: 50
  bg_kmer_cutoff: 10
```

The `columnwise_score_function_name` parameter specifies the function to use for the `columnwise_score_func` argument in the `pairk_conservation_from_json` in `./slim_conservation_scoring/conservation_scores/pairk_conservation.py`. The `bg_cutoff` and `bg_kmer_cutoff` parameters are used as keyword arguments in the function as well.

The current options for the `columnwise_score_function_name` parameter are:
- `shannon_entropy`
- `property_entropy`
These are the same functions as the MSA conservation scores, but they are used to calculate the conservation of the columns in the pairk alignment. To add a different columnwise conservation score to the pipeline, you would import the function into `./slim_conservation_scoring/conservation_scores/__init__.py` and add it to the `ColumnwiseScoreMethods` class as an attribute. You can then use the function in the pipeline by setting the `columnwise_score_function_name` parameter to the attribute name.


## references

- ESM2 (the model used to generate the embeddings): 
    - Z. Lin, H. Akin, R. Rao, B. Hie, Z. Zhu, W. Lu, N. Smetanin, R. Verkuil, O. Kabeli, Y. Shmueli, A. Dos Santos Costa, M. Fazel-Zarandi, T. Sercu, S. Candido, A. Rives, Evolutionary-scale prediction of atomic-level protein structure with a language model. Science 379, 1123–1130 (2023).
- Some of the ESM model sequence encoding functions are adapted from the kibby tool ([link](https://github.com/esbgkannan/kibby)): 
    - W. Yeung, Z. Zhou, S. Li, N. Kannan, Alignment-free estimation of sequence conservation for identifying functional sites using protein sequence embeddings. Brief Bioinform 24 (2023)
- built-in conservation scoring functions are adapted from code released with this study: 
    - J. A. Capra, M. Singh, Predicting functionally important residues from sequence conservation. Bioinformatics 23, 1875–1882 (2007)
- built-in scoring matrix "EDSSMat50" is from this study: 
    - R. Trivedi, H. A. Nagarajaram, Amino acid substitution scoring matrices specific to intrinsically disordered regions in proteins. Sci Rep 9, 16380 (2019)
- built-in "grantham" matrices (including "grantham", "grantham_similarity_norm", and "grantham_similarity_normx100_aligner_compatible") are from or derived from the distance matrix in this study: 
    - R. Grantham, Amino acid difference formula to help explain protein evolution. Science 185, 862–864 (1974).
- other matrices are from biopython:
    - P. J. A. Cock, T. Antao, J. T. Chang, B. A. Chapman, C. J. Cox, A. Dalke, I. Friedberg, T. Hamelryck, F. Kauff, B. Wilczynski, M. J. L. de Hoon, Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics 25, 1422–1423 (2009).




