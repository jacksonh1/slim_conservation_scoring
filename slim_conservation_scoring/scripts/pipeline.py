#!/usr/bin/env python

import argparse
import multiprocessing

from slim_conservation_scoring.pipeline import conservation_pipeline

N_CORES = round(multiprocessing.cpu_count() / 2)


def main_cli():
    parser = argparse.ArgumentParser(
        description="""Run conservation analysis pipeline

provide parameters in a yaml config file (specified with -c)

REQUIRED parameters in the config file:
- database_filekey
- table_file
    - The input table must have a column named 'hit_sequence' if 
        hit_sequence_search_method='search' is used.
    - If hit_sequence_search_method='given_positions', the input table must have
        columns named 'hit start position' and 'hit end position'

Default parameters if not specified in the config file:
output_folder: conservation_analysis
clear_files: True
new_index: True
steps_to_run:
  s1
  s2
  s3
  s4
  s5
  s6
  s7
  s8
  s9
msa_score_methods:
pairk_aln_methods:
pairk_embedding_aln_methods:
hit_sequence_params:
  hit_sequence_search_method: search
  longest_common_subsequence: False
  lcs_min_length: 20
  target_hit_length: 0
idr_params:
  find_idrs: True
  idr_map_file: None
  iupred_cutoff: 0.4
  gap_merge_threshold: 10
  idr_min_length: 8
filter_params:
  min_num_orthos: 20
multilevel_plot_params:
  score_key: aln_property_entropy
  num_bg_scores_cutoff: 20
  score_type: z_score
  strip_gaps: False
aln_slice_params:
  n_flanking_aas: 20
  whole_idr: False
table_annotation_params:
  score_key_for_table: aln_property_entropy
  motif_regex: None
  levels: []
  annotations: ['json_file']
esm_params:
  processes: 4
  threads: 1
  device: cuda
  model_name: esm2_t33_650M_UR50D
pairk_conservation_params:
  kmer_conservation_function_name: pairk_conservation
  columnwise_score_function_name: shannon_entropy
  bg_cutoff: 50
  bg_kmer_cutoff: 10
clean_analysis_files: False
""",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "-c",
        "--config",
        type=str,
        metavar="<file>",
        required=True,
        help="""path to config file in yaml format""",
    )
    parser.add_argument(
        "-n",
        "--n_cores",
        type=int,
        metavar="<int>",
        default=N_CORES,
        help=f"""number of cores to use. Default is {N_CORES}""",
    )
    args = parser.parse_args()
    conservation_pipeline.main(args.config, args.n_cores)


if __name__ == "__main__":
    main_cli()
