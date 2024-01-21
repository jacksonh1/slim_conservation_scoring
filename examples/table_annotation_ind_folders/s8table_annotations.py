# %%
import json
import os
import re
import sys
from pathlib import Path
from typing import Callable

import local_conservation_analysis_pipeline.group_conservation_objects as group_tools
import local_conservation_score_tools.conservation_scoring_tools as cons_tools
import numpy as np
import pandas as pd

# table_file = "./table_original_reindexed.csv"
# levels = ["Eukaryota", "Metazoa", "Vertebrata", "Tetrapoda", "Mammalia"]
# score_key_for_table = "property_entropy"

def get_hit_zscores(lvlo: group_tools.ConserLevel):
    """
    returns a list of the scores and a list of the z scores for the hit (non-gap) positions in query sequence
    """
    hit_slice = slice(lvlo.hit_aln_start, lvlo.hit_aln_end + 1)
    hit_z_scores = lvlo.z_scores[hit_slice]
    hit_scores = lvlo.scores[hit_slice]
    hit_aln_seq = lvlo.query_aln_sequence[hit_slice]
    inds = cons_tools.get_non_gap_indexes(hit_aln_seq)
    return list(np.array(hit_scores)[inds]), list(np.array(hit_z_scores)[inds])

def lvl_annotation_hit_mean_zscore(lvlo: group_tools.ConserLevel):
    scores, z_scores = get_hit_zscores(lvlo)
    return np.mean(z_scores)


def lvl_annotation_aln_slice(lvlo: group_tools.ConserLevel):
    file = Path(lvlo.info_dict["aln_slice_file"]).resolve().relative_to(Path.cwd())
    return rf'=HYPERLINK("{file}")'


def lvl_annotation_conservation_string(lvlo: group_tools.ConserLevel):
    scores, z_scores = get_hit_zscores(lvlo)
    cons_str = cons_tools.conservation_string(
        z_scores, lvlo.hit_aln_sequence.replace('-', ''), z_score_cutoff=0.5
    )
    return cons_str


def addscore(
    json_file: str,
    level: str,
    score_key: str,
    annotation_func: Callable = lvl_annotation_hit_mean_zscore
):
    og = group_tools.ConserGene(json_file)
    if hasattr(og, "critical_error"):
        return
    og.calculate_z_scores(score_key)
    if level not in og.levels_passing_filters:
        return
    if og.level_objects[level].z_score_failure is not None:
        return
    annotation = annotation_func(og.level_objects[level])
    return annotation


def main(table_file, score_key_for_table, levels):
    table_file = table_file.replace(".csv", "_original_reindexed.csv")
    table = pd.read_csv(table_file)
    for level in levels:
        table[f"{level}_{score_key_for_table}_z_score"] = table['json_file'].apply(addscore, args=(level, score_key_for_table, lvl_annotation_hit_mean_zscore))
        table[f"{level}_aln_slice_view"] = table['json_file'].apply(addscore, args=(level, score_key_for_table, lvl_annotation_aln_slice))
        table[f"{level}_cons_string"] = table['json_file'].apply(addscore, args=(level, score_key_for_table, lvl_annotation_conservation_string))
    output_table_file = table_file.replace("_original_reindexed.csv", "_ANALYZED.csv")
    table.to_csv(output_table_file, index=False)
    table.to_excel(output_table_file.replace(".csv", ".xlsx"), index=False)




# %%
'''
def score_obj_2_motif_match_ave_zscore(hits_df, score_o: group_tools.level_score_class):
    ind = hits_df[hits_df["reference_index"] == score_o.reference_index].index[0]
    motif_match = hits_df.loc[ind, "motif_match"]
    full_hit_seq = score_o.hit_sequence_str
    st = full_hit_seq.find(motif_match)
    en = st + len(motif_match) - 1
    hit_z_scores = score_obj_2_nongap_scores(score_o, mask=True)
    return np.mean(hit_z_scores[st : en + 1]


def get_largest_avg_zscore_for_window(score_list, window_size=5):
    """WARNING: this function doesn't account for masked residues (if it's been masked it will count as a 0)"""
    largest_avg_z = np.mean(score_list[0:window_size])
    for i in range(len(score_list) - (window_size - 1)):
        window_scores = score_list[i : i + window_size]
        avg_z = np.mean(window_scores)
        if avg_z > largest_avg_z:
            largest_avg_z = avg_z
    return largest_avg_z


def score_obj_2_largest_avg_zscore_for_window(score_o, window_size=5):
    hit_z_scores = score_obj_2_nongap_scores(score_o)
    largest_avg_z = get_largest_avg_zscore_for_window(
        hit_z_scores, window_size=window_size
    )
    return largest_avg_

def get_MSA_cons_string(score_o):
    hit_z_scores = score_obj_2_nongap_scores(score_o)
    cons_str = cons_tools.conservation_string(
        hit_z_scores, score_o.hit_sequence_str, z_score_cutoff=1
    )
    return cons_str



# %%


def check_json(json_file):
    with open(json_file, "r") as f:
        json_dict = json.load(f)
    if "reference_index" in json_dict:
        return True
    else:
        return False


json_files = Path('./conservation_analysis/').rglob("*.json")
checked_jsons = [i for i in json_files if check_json(i)]

annotation_dict = {}
for i in checked_jsons:
    og = group_tools.ConserGene(i)
    if hasattr(og, "critical_error"):
        continue
    
'''

# %%

# def get_score_for_full_hit(
#     df: pd.DataFrame, 
#     level: str, 
#     score_key: str
# ):
#     for i in df.index:
#         if pd.isna(df.loc[i, "json_file"]):
#             continue
#         og = group_tools.ConserGene(df.loc[i, "json_file"])

#         if level not in og.levels_passing_filters
#             continue
#         z_scores = get_hit_zscores(og.level_objects[level])
#         position_collist = [
#             f"idrbg pos {i} z_score" for i in range(1, len(z_scores) + 1)
#         ]
#         df.loc[i, position_collist] = z_scores
#         df.loc[i, "aln-idrbg-PE-z-score"] = np.mean(z_scores)



# %%

# annotations = {
#     "full_hit": get_score_for_full_hit,
#     "regex match": get_score_for_regex_match_in_hit,
#     "image_file": add_image_file,
#     "alignment_slice": add_alignment_slice,
# }