import json
from pathlib import Path

import local_env_variables.env_variables_and_filepaths as fp
import local_orthoDB_analysis_tools_v2.hit_sequence_tools as hit_tools
import pandas as pd
import pipeline_tools
from pyprojroot import here


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


def main(
    og_info_json: str|Path,
    root: Path,
):
    og_info_json = root / og_info_json
    filter_levels_by_count(i)


