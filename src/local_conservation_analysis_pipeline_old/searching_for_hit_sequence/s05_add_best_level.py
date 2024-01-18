# %%
from math import ceil

import numpy as np

import local_conservation_analysis_pipeline.searching_for_hit_sequence.param_and_log_tools as param_and_log_tools
from local_orthoDB_analysis_tools_v2 import group_tools_v2 as group_tools

# %%
# ==============================================================================
# get the highest ortholog group level with an acceptable number of low gap positions in the hit sequence
# output will be a table column - "orthogroup level"
# ==============================================================================

def enough_usable_residues(score_o, required_usable_fraction=0.5):
    n_residues_required = ceil(len(score_o.hit_sequence_str) * required_usable_fraction)
    hit_slice = slice(
        score_o.hit_alignment_start_position,
        score_o.hit_alignment_end_position
    )
    return np.array(score_o.score_mask[hit_slice]).sum() >= n_residues_required


def find_best_level(og_obj: group_tools.ortholog_group_info, required_usable_fraction=0.5):
    for level in og_obj.levels:
        score_o = og_obj.level_score_objects[level]
        if enough_usable_residues(score_o, required_usable_fraction):
            return level
    return None


def add_best_level_2_json(og_info_file, root, required_usable_fraction=0.5):
    og_obj=group_tools.ortholog_group_info(og_info_file, root=root)
    og_obj.load_all_level_score_objects_with_z_score()
    best_level = find_best_level(og_obj, required_usable_fraction=required_usable_fraction)
    ############ remove or add note for meta json ############
    og_obj.add_item_to_json('best orthogroup level', best_level, save_json=True)
    return og_obj.reference_index, best_level


def main(
    required_usable_fraction,
    meta_info: param_and_log_tools.Metainfo,
    root
):
    og_info_jsons = meta_info.jsons_to_analyze()
    og_info_jsons = [root / i for i in og_info_jsons]
    fail_list = []
    for og_info_json in og_info_jsons:
        ref_ind, best_lvl = add_best_level_2_json(og_info_json, root, required_usable_fraction)
        if best_lvl is None:
            fail_list.append(ref_ind)
    for ref_ind in fail_list:
        meta_info.add_note_for_id(ref_ind, f'no best level found - sequence does not have enough usable residues (>{required_usable_fraction}) in the alignment at any level')
    meta_info.save_meta_dict()
