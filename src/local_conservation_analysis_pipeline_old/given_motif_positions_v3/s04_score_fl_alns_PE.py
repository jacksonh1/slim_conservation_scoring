# %%
import copy
import json
import multiprocessing
from pathlib import Path

import pipeline_tools
from local_orthoDB_analysis_tools_v2 import \
    conservation_scores_2007_capra_singh as cs
from local_orthoDB_analysis_tools_v2 import \
    conservation_scoring_tools as cons_tools
from local_orthoDB_analysis_tools_v2 import group_tools_v2 as group_tools
from pyprojroot import here

# log_file = "./logs/s03.log"
# logger = meta_tools.setup_logging(log_file)


"""
for each protein, 
    load the ortholog group info json
    for each level
        load the alignment info json
        load the alignment
        gap_mask, score_mask = create the mask using: query aligned sequence, gap_frac_cutoff
        calculate the scores using: query aligned sequence, alignment, matrix_df
        save the scores to a json
        add the score_json file name to the alignment info
"""


def level_info_2_mask(level_info: group_tools.level_info, gap_frac_cutoff=0.2):
    gap_mask, score_mask = cons_tools.make_score_mask(
        level_info.alignment_orthologs_ldo_clustered,
        level_info.query_aligned_sequence_str,
        gap_frac_cutoff=gap_frac_cutoff,
    )
    gap_mask = [bool(i) for i in gap_mask]
    score_mask = [bool(i) for i in score_mask]
    return gap_mask, score_mask


def level_info_2_scores(level_info: group_tools.level_info):
    aln = copy.deepcopy(level_info.alignment_orthologs_ldo_clustered)
    scores = []
    for i in range(aln.get_alignment_length()):
        col = aln[:, i]
        scores.append(cs.property_entropy(col))
    return scores


def level_info_2_score_dict(
    level_info: group_tools.level_info, gap_frac_cutoff=0.2
):
    score_dict = {}
    score_dict["gap_mask"], score_dict["score_mask"] = level_info_2_mask(
        level_info, gap_frac_cutoff=gap_frac_cutoff
    )
    score_dict["scores"] = level_info_2_scores(level_info)
    return score_dict


def level_info_2_score_json(
    level_info: group_tools.level_info, gap_frac_cutoff=0.2
):
    score_dict = level_info_2_score_dict(
        level_info, gap_frac_cutoff=gap_frac_cutoff
    )
    score_json = (
        level_info.analysis_folder
        / f"conservation_scores_property_entropy_{level_info.level}.json"
    )
    with open(score_json, "w") as f:
        json.dump(score_dict, f, indent=4)
    return score_json


def score_levels_in_og_obj(
    og_obj: group_tools.ortholog_group_info, gap_frac_cutoff, root
):
    for level in og_obj.levels:
        level_info = og_obj.level_info_objects[level]
        if og_obj.info_dict["num_sequences_per_level"][level] < 10:
            # og_obj.add_item_to_level_info(
                # "property_entropy_score_json",
                # "fewer than 10 sequences",
                # level,
            # )
            continue
        score_json = level_info_2_score_json(
            level_info, gap_frac_cutoff=gap_frac_cutoff
        )
        og_obj.add_item_to_level_info(
            "property_entropy_score_json",
            str(score_json.resolve().relative_to(root)),
            level,
        )


def driver(og_info_json, gap_frac_cutoff, root):
    # ic(og_info_json)
    og_obj = group_tools.ortholog_group_info(og_info_json, root=root)
    og_obj.load_all_level_info_objects()
    score_levels_in_og_obj(og_obj, gap_frac_cutoff=gap_frac_cutoff, root=root)


def main(
    gap_frac_cutoff,
    og_info_json: str|Path,
    root: Path,
):
    og_info_json = root / og_info_jsons
    driver(og_info_json, gap_frac_cutoff, root)
