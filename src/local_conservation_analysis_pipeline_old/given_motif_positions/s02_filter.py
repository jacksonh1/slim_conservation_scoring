import json
from pathlib import Path

import pandas as pd
from pyprojroot import here

import local_conservation_analysis_pipeline.given_motif_positions.param_and_log_tools as param_and_log_tools
import local_env_variables.env_variables_and_filepaths as fp
import local_orthoDB_analysis_tools_v2.hit_sequence_tools as hit_tools


def levels_above_cutoff(level_counts_dict, count_cutoff):
    return [i for i in level_counts_dict.keys() if level_counts_dict[i] >= count_cutoff]

def filter_levels_by_count(og_info_json):
    with open(og_info_json, 'r') as f:
        og_info = json.load(f)
    level_counts_dict = og_info["num_sequences_per_level"]
    levels = og_info["levels"]
    levels_with_enough_seqs = levels_above_cutoff(level_counts_dict, 10)
    # using nested list to preserve order
    levels = [i for i in levels if i in levels_with_enough_seqs]
    og_info["levels"] = levels
    with open(og_info_json, 'w') as f:
        json.dump(og_info, f, indent=4)

def remove_failed_hits_from_meta_dict(
    og_info_json,
    meta_info: param_and_log_tools.Metainfo
):
    with open(og_info_json, 'r') as f:
        og_info = json.load(f)
    if not og_info['found hit']:
        meta_info.add_id_to_fail_dict(og_info["reference_index"], 'hit sequence not found')
        return
    if len(og_info['levels']) == 0:
        meta_info.add_id_to_fail_dict(og_info["reference_index"], f'no levels with >=10 orthologs')
        return
    if not og_info['hit_in_idr']:
        meta_info.add_id_to_fail_dict(og_info["reference_index"], 'hit sequence found but it is not in an IDR')
        return

def main(
    meta_info: param_and_log_tools.Metainfo,
    root: Path,
):
    og_info_jsons = meta_info.jsons_to_analyze()
    og_info_jsons = [root / i for i in og_info_jsons]
    for i in og_info_jsons:
        filter_levels_by_count(i)
        remove_failed_hits_from_meta_dict(i, meta_info)
    meta_info.save_meta_dict()


