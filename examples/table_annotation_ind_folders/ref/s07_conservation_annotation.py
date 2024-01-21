# %%
"""
add conservation string
add sequence distance z-score?

"""
import json
import traceback

import meta_tools
import numpy as np
import pandas as pd
from local_orthoDB_analysis_tools_v2 import \
    conservation_scoring_tools as cons_tools
from local_orthoDB_analysis_tools_v2 import group_tools_v2 as group_tools
from pyprojroot import here

# %load_ext autoreload
# %autoreload 2
HITS_DF = pd.read_csv(meta_tools.HITS_FILE.replace(".csv", "_original_reindexed.csv"))
LEVEL_TO_ANALYZE = meta_tools.params["LEVEL_TO_USE_FOR_SCORES"]
SCORE_KEY = meta_tools.params["SCORE_TYPE_FOR_TABLE_KEY"]

import local_env_variables.env_variables_and_filepaths as fp

ROOT = fp.conservation_analysis_root

# %%
# ==============================================================================
# get some quantities to annotate the table with
# full MSA conservation scores:
#   - level looked at
#   - conservation string
#   - average z-score for best 5 residue window
#   - average z-score for best 10 residue window
#   - average z-score across full window
#   - # points in z-score background?
# save all of this to a json file
# ==============================================================================
# %%


# %%
def score_obj_2_motif_match_ave_zscore(hits_df, score_o: group_tools.level_score_class):
    ind = hits_df[hits_df["reference_index"] == score_o.reference_index].index[0]
    motif_match = hits_df.loc[ind, "motif_match"]
    full_hit_seq = score_o.hit_sequence_str
    st = full_hit_seq.find(motif_match)
    en = st + len(motif_match) - 1
    hit_z_scores = score_obj_2_nongap_scores(score_o, mask=True)
    return np.mean(hit_z_scores[st : en + 1])


# %%


def mask_score_list(score_list, score_mask):
    score_list = np.array(score_list)
    score_list[[not i for i in score_mask]] = 0
    return list(score_list)


def score_obj_2_nongap_scores(score_o: group_tools.level_score_class, mask=True):
    if mask:
        z_score_list = mask_score_list(
            score_o.z_score_dict["z_scores"], score_o.score_mask
        )
    else:
        z_score_list = score_o.z_score_dict["z_scores"]
    hit_aln_slice = slice(
        score_o.hit_alignment_start_position, score_o.hit_alignment_end_position + 1
    )
    nongap_positions = cons_tools.get_non_gap_indexes(
        score_o.query_aligned_sequence_str[hit_aln_slice]
    )
    hit_z_scores = np.array(z_score_list[hit_aln_slice])[nongap_positions]
    return list(hit_z_scores)


def get_MSA_cons_string(score_o):
    hit_z_scores = score_obj_2_nongap_scores(score_o)
    cons_str = cons_tools.conservation_string(
        hit_z_scores, score_o.hit_sequence_str, z_score_cutoff=1
    )
    return cons_str


# write a function to get the largest average z-score for a window of residues
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
    return largest_avg_z


# %%


def append_dict(score_o: group_tools.level_score_class, annot_dict={}):
    if score_o.z_score_dict is None:
        return annot_dict
    annot_dict[score_o.reference_index] = {}
    annot_dict[score_o.reference_index]["full MSA - OG level"] = score_o.level
    annot_dict[score_o.reference_index]["full MSA - cons str"] = get_MSA_cons_string(
        score_o
    )
    annot_dict[score_o.reference_index][
        "full MSA - average score of best-scoring 5-residue window"
    ] = score_obj_2_largest_avg_zscore_for_window(score_o, window_size=5)
    annot_dict[score_o.reference_index][
        "full MSA - average score of best-scoring 10-residue window"
    ] = score_obj_2_largest_avg_zscore_for_window(score_o, window_size=10)
    annot_dict[score_o.reference_index][
        "full MSA - average score over full window"
    ] = score_obj_2_largest_avg_zscore_for_window(
        score_o, window_size=len(score_o.hit_sequence_str)
    )
    annot_dict[score_o.reference_index][
        "full MSA - # points in z-score background"
    ] = len(score_o.z_score_dict["bg_scores"])
    annot_dict[score_o.reference_index][
        "full MSA - z-score list"
    ] = score_obj_2_nongap_scores(score_o, mask=False)
    # annot_dict[score_o.reference_index]['full MSA - motif_match_score'] = score_obj_2_motif_match_ave_zscore(HITS_DF, score_o)
    return annot_dict


def driver(
    ortholog_group_info_json,
    level="best",
    annot_dict={},
    local_sequences=False,
    score_key="asym_valday_score_json",
):
    print(ortholog_group_info_json)
    og_obj = group_tools.ortholog_group_info(ortholog_group_info_json, root=ROOT)
    assert level in og_obj.levels or level == "best"
    if level == "best":
        if og_obj.best_level is None:
            return annot_dict
        else:
            level = og_obj.best_level
    if not og_obj.hit_in_idr:
        return annot_dict
    try:
        score_o = og_obj.load_level_score_object_w_z_score(
            level,
            local_sequences=local_sequences,
            score_key=score_key,
        )
        # og_obj.load_all_level_score_objects_with_z_score(local_sequences=True)
    except Exception as e:
        traceback.print_exc()
        print(e)
        return annot_dict
    # score_o = og_obj.level_score_objects[level]
    annot_dict = append_dict(score_o, annot_dict=annot_dict)
    return annot_dict


def main():
    og_info_jsons = meta_tools.jsons_to_analyze()
    og_info_jsons = [ROOT / i for i in og_info_jsons]
    annotation_dict = {}
    for og_info_json in og_info_jsons:
        annotation_dict = driver(
            og_info_json, LEVEL_TO_ANALYZE, annotation_dict, score_key=SCORE_KEY
        )
    with open("./s07_full_msa_annotation.json", "w") as f:
        json.dump(annotation_dict, f, indent=4)


if __name__ == "__main__":
    main()
