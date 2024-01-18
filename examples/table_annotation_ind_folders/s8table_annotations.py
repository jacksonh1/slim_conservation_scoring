import json
import os
import re
import sys

import local_conservation_analysis_pipeline.group_conservation_objects as group_tools
import local_conservation_score_tools.conservation_scoring_tools as cons_tools
import numpy as np
import pandas as pd

table_file = './table_original_reindexed.csv'
table = pd.read_csv(table_file)


def get_hit_zscores(score_o: group_tools.ConserLevel, score_key: str):
    hit_slice = slice(score_o.hit_alignment_start_position, score_o.hit_alignment_end_position+1)
    hit_scores = score_o.z_score_dict['z_scores'][hit_slice]
    hit_aln_seq = score_o.query_aligned_sequence_str[hit_slice]
    inds = cons_tools.get_non_gap_indexes(hit_aln_seq)
    return list(np.array(hit_scores)[inds])


def get_score_for_full_hit():
    for i in df.index:
        og = group_tools.ortholog_group_info(ROOT / df.loc[i, 'og_info_json'])
        level = 'Vertebrata'
        try:
            score_o = og.load_level_score_object_w_z_score(level, local_sequences=True, score_key='property_entropy_score_json')
        except KeyError as ke:
            print(f"{df.loc[i, 'reference_index']} - {ke}")
            continue
        if score_o.z_score_dict is None:
            print(f"no z scores for {df.loc[i, 'reference_index']} - {score_o.z_score_failure_reason}")
            continue
        z_scores = get_hit_zscores(score_o)
        position_collist = [f'idrbg pos {i} z_score' for i in range(1,len(z_scores)+1)]
        df.loc[i, position_collist] = z_scores
        df.loc[i, 'aln-idrbg-PE-z-score'] = np.mean(z_scores)

annotations = {
    "full_hit": get_score_for_full_hit,
    "regex match": get_score_for_regex_match_in_hit,
    "image_file": add_image_file,
    "alignment_slice": add_alignment_slice,
}




