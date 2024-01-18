"""
for when you only have the uniprot id and the hit sequence
What is the goal of this script?
- create an index representing each row of the hits table
- create folders for the ortholog analysis with a folder for each table entry (index)
- find the orthoDB id for each uniprot id
- Access the orthoDB json files for each level (files from a previous orthogroup construction run `odb_og_run_folder`)
- create a json file for storing information about the ortholog groups for each table entry
- add some info to that json file
- log table entries that failed to be processed (most likely because no ortholog 
group folder was found or the uniprot id couldn't be mapped to an orthoDB id)
"""

import json
from pathlib import Path

import pandas as pd

import local_conservation_analysis_pipeline.searching_for_hit_sequence.param_and_log_tools as param_and_log_tools
import local_env_variables.env_variables_and_filepaths as fp
from local_orthoDB_tools import uniprotid_search


def _get_number_seqs_in_level(odb_json_file):
    with open(odb_json_file, "r") as f:
        d = json.load(f)
    return d["num_sequences_LDO_cdhit"]


def get_number_seqs_in_all_levels(level_jsons_dict):
    levels_seq_count_dict = {
        i: _get_number_seqs_in_level(level_jsons_dict[i])
        for i in level_jsons_dict.keys()
    }
    return levels_seq_count_dict


def _get_og_id_for_level(odb_json_file):
    with open(odb_json_file, "r") as f:
        d = json.load(f)
    return d["selected_query_ogid"]


def get_og_ids_in_all_levels(level_jsons_dict):
    levels_og_ids_dict = {
        i: _get_og_id_for_level(level_jsons_dict[i]) for i in level_jsons_dict.keys()
    }
    return levels_og_ids_dict


def import_and_reindex_hits_df(hits_file):
    hits_df = pd.read_csv(hits_file)
    hits_df = hits_df.reset_index(names="reference_index")
    hits_df.to_csv(hits_file.replace(".csv", "_original_reindexed.csv"), index=False)
    return hits_df


def get_geneid_folder_dict(odb_main_folder):
    odb_main_folder = Path(odb_main_folder)
    return {x.name: x for x in odb_main_folder.glob("*/") if x.is_dir()}


def get_levels_json_dict(odb_folder):
    odb_folder = Path(odb_folder)
    levels_json_dict = {str(i.parent.name): str(i) for i in odb_folder.glob("*/*.json")}
    return levels_json_dict


def get_sorted_levels_list(levels_json_dict):
    df = fp.load_data_levels_df()
    df = df[df["level name"].isin(levels_json_dict.keys())]
    df = df.astype({"total non-redundant count of species underneath": "int"})
    df = df.sort_values(
        "total non-redundant count of species underneath", ascending=False
    )
    return list(df["level name"])


def ortholog_analysis_setup(
    hits_df, odb_main_folder, og_analysis_output_folder, root, logger
):
    odb_main_folder = Path(odb_main_folder)
    geneid_folders_dict = get_geneid_folder_dict(odb_main_folder)

    references_failed = {}
    references_passed = {}
    for i in hits_df.index:
        id_o = hits_df.loc[i, "UniprotID"]
        ref_ind = int(hits_df.loc[i, "reference_index"])
        try:
            logger.info(f"searching for {id_o} in orthoDB")
            geneid = uniprotid_search.uniprotid_2_geneid_human(id_o)
        except ValueError as ve:
            logger.warning(f"no gene id found for {id_o} - {ve}")
            logger.warning(ve, exc_info=True)
            # hits_df.loc[i, 'odb_geneid'] = False
            references_failed[ref_ind] = "no odb gene id found for uniprot id"
            continue
        try:
            odb_folder = str(geneid_folders_dict[geneid])
        except KeyError as ke:
            logger.warning(
                f"no ortholog folder found for {id_o} (gene_id: {geneid}) - {ke}"
            )
            # hits_df.loc[i, 'odb_geneid'] = False
            references_failed[ref_ind] = "no ortholog folder found for gene id"
            continue

        hit_seq = hits_df.loc[i, "hit sequence"]
        if "name" in hits_df.columns:
            name = hits_df.loc[i, "name"]
        else:
            name = id_o
        if "regex" in hits_df.columns:
            regex = hits_df.loc[i, "regex"]
        else:
            regex = None
        ref_ind = int(hits_df.loc[i, "reference_index"])
        level_jsons = get_levels_json_dict(odb_folder)
        levels_list = get_sorted_levels_list(level_jsons)
        level_counts = get_number_seqs_in_all_levels(level_jsons)
        level_og_ids = get_og_ids_in_all_levels(level_jsons)

        info_dict = {
            "UniprotID": id_o,
            "query_gene_id": geneid,
            "odb_folder": odb_folder,
            "orig_hit_sequence": hit_seq,
            "name": name,
            "reference_index": ref_ind,
            "json_files_per_level": level_jsons,
            "num_sequences_per_level": level_counts,
            "levels_og_ids_dict": level_og_ids,
            # 'level_files': {},
            "levels": levels_list,
            "regex": regex,
        }

        # for level in levels_list:
        # info_dict['level_files'][level] = {}
        # info_dict['level_files'][level]['odb_json'] = level_jsons[level]

        hit_output_folder = og_analysis_output_folder / f"{ref_ind}-{id_o}_{geneid}"
        hit_output_folder.mkdir(exist_ok=True, parents=True)
        group_info_json = hit_output_folder / "ortholog_group_info.json"
        with open(group_info_json, "w") as f:
            json.dump(info_dict, f, indent=4)

        references_passed[ref_ind] = str(group_info_json.resolve().relative_to(root))
    return hits_df, references_failed, references_passed


def main(
    hits_file,
    og_analysis_output_folder,
    odb_og_run_folder,
    logger,
    meta_info: param_and_log_tools.Metainfo,
    root,
):
    meta_info.clear_meta_json()
    og_analysis_output_folder = Path(og_analysis_output_folder)
    og_analysis_output_folder.mkdir(exist_ok=True, parents=True)
    og_analysis_output_folder_rel_to_root = (
        og_analysis_output_folder.resolve().relative_to(root)
    )
    hits_df = import_and_reindex_hits_df(hits_file)
    (
        hits_df,
        references_failed,
        references_passed,
    ) = ortholog_analysis_setup(
        hits_df, odb_og_run_folder, og_analysis_output_folder, root, logger
    )
    meta_info.META_PARAMS["og_analysis_output_folder"] = str(
        og_analysis_output_folder_rel_to_root
    )
    meta_info.META_PARAMS["odb_main_folder"] = str(odb_og_run_folder)
    meta_info.META_PARAMS["ref_id_dict"] = references_passed
    meta_info.META_PARAMS["failed_id_dict"] = references_failed
    meta_info.META_PARAMS["notes"] = {}
    meta_info.save_meta_dict()
