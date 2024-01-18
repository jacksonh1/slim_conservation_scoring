# %%
import json
import multiprocessing

from Bio import AlignIO, Seq, SeqIO
from pyprojroot import here

import local_conservation_analysis_pipeline.given_motif_positions.param_and_log_tools as param_and_log_tools
import local_env_variables.matrices as submats
from local_orthoDB_analysis_tools_v2 import \
    conservation_scoring_tools as cons_tools
from local_orthoDB_analysis_tools_v2 import group_tools_v2 as group_tools
from pathlib import Path

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


def level_info_2_scores(level_info: group_tools.level_info, matrix_df):
    query_record = SeqIO.SeqRecord(seq=Seq.Seq(level_info.query_aligned_sequence_str), id=level_info.query_sequence_id_str, description='')
    scores = cons_tools.asymmetric_valdar_score_df_unweighted(
        query_record,
        level_info.alignment_orthologs_ldo_clustered,
        matrix_df,
    )
    return scores


def level_info_2_score_dict(
    level_info: group_tools.level_info, matrix_df, gap_frac_cutoff=0.2
):
    score_dict = {}
    score_dict["gap_mask"], score_dict["score_mask"] = level_info_2_mask(
        level_info, gap_frac_cutoff=gap_frac_cutoff
    )
    score_dict["scores"] = level_info_2_scores(level_info, matrix_df)
    return score_dict


def level_info_2_score_json(
    level_info: group_tools.level_info,
    matrix_df,
    matrix_name,
    output_folder,
    gap_frac_cutoff,
):
    score_json = (
        Path(output_folder)
        / f"{level_info.query_gene_id}_conservation_scores_asym_valday_{level_info.level}_{matrix_name}.json"
    )
    if score_json.exists():
        print(f"score_json {score_json} already exists")
        return score_json
    score_dict = level_info_2_score_dict(
        level_info, matrix_df, gap_frac_cutoff=gap_frac_cutoff
    )
    with open(score_json, "w") as f:
        json.dump(score_dict, f, indent=4)
    return score_json


def score_levels_in_og_obj(
    og_obj: group_tools.ortholog_group_info,
    matrix_name,
    gap_frac_cutoff,
    output_folder,
    root
):
    matrix_df = submats.load_precomputed_matrix_df(matrix_name)
    for level in og_obj.levels:
        level_info = og_obj.level_info_objects[level]
        if og_obj.info_dict["num_sequences_per_level"][level] < 10:
            continue
        score_json = level_info_2_score_json(
            level_info,
            matrix_df,
            matrix_name,
            output_folder=output_folder,
            gap_frac_cutoff=gap_frac_cutoff
        )
        og_obj.add_item_to_level_info(
            "asym_valday_score_json",
            str(score_json.resolve().relative_to(root)),
            level,
        )

def driver(og_info_json, matrix_name, gap_frac_cutoff, output_folder, root):
    og_obj = group_tools.ortholog_group_info(og_info_json, root=root)
    og_obj.load_all_level_info_objects()
    score_levels_in_og_obj(og_obj, matrix_name, gap_frac_cutoff, output_folder, root)

def main(
    matrix_name,
    gap_frac_cutoff,
    output_folder,
    root,
    meta_info: param_and_log_tools.Metainfo,
    multiprocess=True
):
    og_info_jsons = meta_info.jsons_to_analyze()
    og_info_jsons = [root / i for i in og_info_jsons]
    output_folder = Path(output_folder)
    output_folder.mkdir(exist_ok=True, parents=True)
    if multiprocess:
        p = multiprocessing.Pool(multiprocessing.cpu_count() - 1)
        f_args = [(i, matrix_name, gap_frac_cutoff, output_folder, root) for i in og_info_jsons]
        p.starmap(driver, f_args)
        # might be better to use imap_unordered or something that will close each process as it finishes
        # p.map(driver, og_info_jsons)
        p.close()
        p.join()
    else:
        # driver([i for i in og_info_jsons if '9606_0:0014fb' in i][0], PARAMS)
        for i in og_info_jsons:
            driver(i, matrix_name, gap_frac_cutoff, output_folder, root)

