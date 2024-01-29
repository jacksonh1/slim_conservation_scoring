import json
from pathlib import Path

import numpy as np
import pandas as pd
from openpyxl import Workbook, load_workbook

import local_conservation_analysis_pipeline.group_conservation_objects as group_tools
import local_conservation_score_tools.score_tools as cons_tools
import local_seqtools.general_utils as tools


def check_json(json_file):
    with open(json_file, "r") as f:
        json_dict = json.load(f)
    if "reference_index" in json_dict:
        return True
    else:
        return False


def get_image_file(og: group_tools.ConserGene, image_score_key):
    if f"multilevel_plot_file-{image_score_key}" in og.info_dict:
        file = Path(og.info_dict[f"multilevel_plot_file-{image_score_key}"]).resolve()
        # file = Path(og.info_dict[f"multilevel_plot_file-{image_score_key}"]).resolve().relative_to(Path.cwd())
        return rf'=HYPERLINK("./{file}")'


def find_motif_regex(og: group_tools.ConserGene, regex):
    if regex is None:
        return None, None, None
    hit_sequence = og.hit_sequence
    matches = list(tools.get_regex_matches(regex, hit_sequence))
    if len(matches) == 0:
        return None, None, None
    if len(matches) > 1:
        return None, None, None
    m = matches[0]
    matchseq = m[0]
    matchst = m[1]
    matchend = m[2]
    return matchseq, matchst, matchend


def get_hit_aln_scores(lvlo: group_tools.LevelAlnScore):
    """
    returns a list of the scores and a list of the z scores for the hit (non-gap) positions in query sequence
    """
    hit_slice = slice(lvlo.hit_aln_start, lvlo.hit_aln_end + 1)
    hit_z_scores = lvlo.z_scores[hit_slice]
    hit_scores = lvlo.scores[hit_slice]
    hit_aln_seq = lvlo.query_aln_sequence[hit_slice]
    nongap_inds = tools.get_non_gap_indexes(hit_aln_seq)
    return list(np.array(hit_scores)[nongap_inds]), list(
        np.array(hit_z_scores)[nongap_inds]
    )


def lvl_annotation_conservation_string(lvlo: group_tools.LevelAlnScore):
    scores, z_scores = get_hit_aln_scores(lvlo)
    cons_str = cons_tools.conservation_string(
        z_scores, lvlo.hit_aln_sequence.replace('-', ''), z_score_cutoff=0.5
    )
    return cons_str


def lvl_annotation_aln_slice(lvlo: group_tools.LevelAlnScore):
    file = Path(lvlo.info_dict["aln_slice_file"]).resolve().relative_to(Path.cwd())
    return rf'=HYPERLINK("{file}")'


def get_largest_avg_zscore_for_window(score_list, window_size=5):
    largest_avg_z = np.mean(score_list[0:window_size])
    for i in range(len(score_list) - (window_size - 1)):
        window_scores = score_list[i : i + window_size]
        avg_z = np.mean(window_scores)
        if avg_z > largest_avg_z:
            largest_avg_z = avg_z
    return largest_avg_z


def add_gene_annotations_2_dict(
    json_file, annotation_dict, image_score_key, table_annotation_score_key, regex=None
):
    og = group_tools.ConserGene(json_file)
    ref_ind = og.reference_index
    if ref_ind in annotation_dict:
        raise ValueError(f"reference index {ref_ind} already in annotation_dict")
    d = {}
    d["json_file"] = str(json_file.resolve())
    if "critical_error" in og.info_dict:
        d["critical_error"] = og.info_dict["critical_error"]
        annotation_dict[ref_ind] = d
        return annotation_dict
    else:
        d["critical_error"] = None
    d["hit_start_position"] = og.hit_start_position
    # d['hit_end_position'] = og.hit_end_position
    d["multi_level_plot"] = get_image_file(og, image_score_key)
    match = find_motif_regex(og, regex)
    d["regex"] = regex
    d["regex_match"] = match[0]
    d["regex_match_stpos_in_hit"] = match[1]
    d["regex_match_endpos_in_hit"] = match[2]
    # alignment score annotations
    og.load_aln_scores(score_key=table_annotation_score_key)
    d['level_annotations'] = {}
    for level, lvlo in og.aln_score_objects.items():
        if level in d['level_annotations']:
            raise ValueError(f"level {level} already in d['level_annotations']")
        d['level_annotations'][level] = {}
        dlvl = d['level_annotations'][level]
        hit_scores, hit_z_scores = get_hit_aln_scores(lvlo)
        dlvl["hit_z_scores"] = hit_z_scores
        dlvl["hit_mean_zscore"] = np.mean(hit_z_scores)
        dlvl["hit_scores"] = hit_scores
        dlvl["hit_mean_score"] = np.mean(hit_scores)
        dlvl["aln_slice_file"] = lvl_annotation_aln_slice(lvlo)
        dlvl["conservation_string"] = lvl_annotation_conservation_string(lvlo)
        if len(hit_z_scores) >= 5:
            dlvl["best mean z-score over 5 residue window"] = get_largest_avg_zscore_for_window(
                hit_z_scores, window_size=5
            )
        if len(hit_z_scores) >= 10:
            dlvl["best mean z-score over 10 residue window"] = get_largest_avg_zscore_for_window(
                hit_z_scores, window_size=10
            )
        if match[1] is not None:
            regex_scores = hit_scores[match[1] : match[2] + 1]
            dlvl["regex_match_scores"] = regex_scores
            dlvl["regex_match_mean_score"] = np.mean(regex_scores)
            regex_z_scores = hit_z_scores[match[1] : match[2] + 1]
            dlvl["regex_match_zscores"] = regex_z_scores
            dlvl["regex_match_mean_zscore"] = np.mean(regex_z_scores)
    annotation_dict[ref_ind] = d
    return annotation_dict


def main(search_dir, image_score_key, table_annotation_score_key, regex=None):
    json_files = Path(search_dir).rglob("*.json")
    checked_jsons = [i for i in json_files if check_json(i)]
    annotations = {}
    for json_file in checked_jsons:
        annotations = add_gene_annotations_2_dict(
            json_file,
            annotations,
            image_score_key,
            table_annotation_score_key,
            regex=regex,
        )
    # print(checked_jsons)
    # print(annotations)
    # return annotations
    if Path("annotations.json").exists():
        # remove old annotations file
        Path("annotations.json").unlink()
    with open("annotations.json", "w") as f:
        json.dump(annotations, f, indent=4)