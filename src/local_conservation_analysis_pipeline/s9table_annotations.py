
import json
import os
import re
import sys
from functools import partial
from pathlib import Path
from typing import Callable

import numpy as np
import pandas as pd

import local_conservation_analysis_pipeline.group_conservation_objects as group_tools
import local_conservation_score_tools.score_tools as cons_tools
import local_seqtools.general_utils as tools

# table_file = "./table_original_reindexed.csv"
# levels = ["Eukaryota", "Metazoa", "Vertebrata", "Tetrapoda", "Mammalia"]
# score_key_for_table = "property_entropy"

def get_hit_zscores(lvlo: group_tools.LevelAlnScore):
    """
    returns a list of the scores and a list of the z scores for the hit (non-gap) positions in query sequence
    """
    hit_slice = slice(lvlo.hit_aln_start, lvlo.hit_aln_end + 1)
    hit_z_scores = lvlo.z_scores[hit_slice]
    hit_scores = lvlo.scores[hit_slice]
    hit_aln_seq = lvlo.query_aln_sequence[hit_slice]
    inds = tools.get_non_gap_indexes(hit_aln_seq)
    return list(np.array(hit_scores)[inds]), list(np.array(hit_z_scores)[inds])


def lvl_annotation_hit_mean_zscore(lvlo: group_tools.LevelAlnScore):
    scores, z_scores = get_hit_zscores(lvlo)
    return np.mean(z_scores)


def lvl_annotation_hit_mean_score(lvlo: group_tools.LevelAlnScore):
    scores, z_scores = get_hit_zscores(lvlo)
    return np.mean(scores)


def lvl_annotation_aln_slice(lvlo: group_tools.LevelAlnScore):
    file = Path(lvlo.info_dict["aln_slice_file"]).resolve().relative_to(Path.cwd())
    return rf'=HYPERLINK("{file}")'


def lvl_annotation_conservation_string(lvlo: group_tools.LevelAlnScore):
    scores, z_scores = get_hit_zscores(lvlo)
    cons_str = cons_tools.conservation_string(
        z_scores, lvlo.hit_aln_sequence.replace('-', ''), z_score_cutoff=0.5
    )
    return cons_str


def lvl_annotation_regex_match(lvlo: group_tools.LevelAlnScore, regex: str):
    hit_slice = slice(lvlo.hit_aln_start, lvlo.hit_aln_end + 1)
    hit_seq = lvlo.query_aln_sequence[hit_slice].replace('-', '')
    scores, z_scores = get_hit_zscores(lvlo)
    matches = list(tools.get_regex_matches(regex, hit_seq))
    if len(matches) == 0:
        return
    if len(matches) > 1:
        return
    m = matches[0]
    matchst = m[1]
    matchen = m[2]
    match_z_scores = z_scores[matchst:matchen + 1]
    return m[0], np.mean(match_z_scores)
'''
Might be better to at a column to the table with the regex match and it's relative position in the hit sequence. Then annotate the table with the z-score of the regex match. Could use jch_alignment to slice up the alignment and convert positions
'''

def addscore(
    json_file: str,
    level: str,
    score_key: str,
    annotation_func: Callable = lvl_annotation_hit_mean_zscore
):
    og = group_tools.ConserGene(json_file)
    if hasattr(og, "critical_error"):
        return
    og.load_aln_scores(score_key)
    if level not in og.levels_passing_filters:
        return
    if og.aln_score_objects[level].z_score_failure is not None:
        return
    annotation = annotation_func(og.aln_score_objects[level])
    return annotation


def main(table_file, score_key_for_table, levels, regex=None):
    table_file = table_file.replace(".csv", "_original_reindexed.csv")
    table = pd.read_csv(table_file)
    for level in levels:
        table[f"{level}_{score_key_for_table}_z_score"] = table['json_file'].apply(addscore, args=(level, score_key_for_table, lvl_annotation_hit_mean_zscore))
        table[f"{level}_aln_slice_view"] = table['json_file'].apply(addscore, args=(level, score_key_for_table, lvl_annotation_aln_slice))
        table[f"{level}_cons_string"] = table['json_file'].apply(addscore, args=(level, score_key_for_table, lvl_annotation_conservation_string))
    if regex is not None:
        lvl_annotation_regex_match_func = partial(lvl_annotation_regex_match, regex=regex)
        for level in levels:
            # add results as 2 columns: regex_match, regex_match_z_score
            table[f"{level}_regex_match"] = table['json_file'].apply(addscore, args=(level, score_key_for_table, lvl_annotation_regex_match_func))
            table[f"{level}_regex_match_z_score"] = table[f"{level}_regex_match"].apply(lambda x: x[1] if x is not None else None)
            table[f"{level}_regex_match"] = table[f"{level}_regex_match"].apply(lambda x: x[0] if x is not None else None)
    output_table_file = table_file.replace("_original_reindexed.csv", "_ANALYZED.csv")
    table.to_csv(output_table_file, index=False)
    table.to_excel(output_table_file.replace(".csv", ".xlsx"), index=False)



