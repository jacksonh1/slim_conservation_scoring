"""
Write a separate script to find and annotate the table with the odbid
for when you already have the orthodb id and the hit sequence positions in the orthodb sequence
goal of this script:
- create an index representing each row of the hits table
- create folders for the ortholog analysis with a folder for each table entry (index)
- Access the orthoDB json files for each level (files from a previous orthogroup construction run `odb_og_run_folder`)
- create a json file for storing information about the ortholog groups for each table entry
- add some info to that json file
- log table entries that failed to be processed (most likely because no ortholog 
group folder was found)
"""
import json
from pathlib import Path

# import meta_tools
import pandas as pd
from pyprojroot import here

import local_conservation_analysis_pipeline.given_motif_positions.param_and_log_tools as param_and_log_tools
import local_env_variables.env_variables_and_filepaths as fp
import local_orthoDB_analysis_tools_v2.hit_sequence_tools as hit_tools

# logger = meta_tools.setup_logging("./logs/s01.log")
# HITS_FILE = meta_tools.HITS_FILE
# ODB_MAIN_FOLDER = meta_tools.ODB_MAIN_FOLDER
# OG_ANALYSIS_OUTPUT_FOLDER = meta_tools.OG_ANALYSIS_OUTPUT_FOLDER
# META_FILE = meta_tools.META_FILE
# ROOT = fp.conservation_analysis_root


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


def get_sorted_levels_list(levels_json_dict):
    df = fp.load_data_levels_df()
    df = df[df["level name"].isin(levels_json_dict.keys())]
    df = df.astype({"total non-redundant count of species underneath": "int"})
    df = df.sort_values(
        "total non-redundant count of species underneath", ascending=False
    )
    return list(df["level name"])


def create_level_info_dict(
    levels, level_jsons_dict, hit_start_position, hit_end_position, **extra_info
):
    level_info_dict = {}
    for level in levels:
        try:
            lvl_info_obj = hit_tools.hit_sequence_class_given_hit_positions(
                level_jsons_dict[level], hit_start_position, hit_end_position
            )
        except IndexError as ie:
            raise
        level_info_dict[level] = lvl_info_obj.get_attributes_as_dict()
        level_info_dict[level].update(extra_info)
        # breakpoint()
    return level_info_dict


def add_lvlinfo_2_info_dict(info_dict, level_info_dict):
    info_dict["query_sequence_id_str"] = level_info_dict["query_sequence_id_str"]
    info_dict["query_sequence_str"] = level_info_dict["query_sequence_str"]
    info_dict["hit_sequence"] = level_info_dict["hit_sequence"]
    info_dict["hit_start_position"] = level_info_dict["hit_start_position"]
    info_dict["hit_end_position"] = level_info_dict["hit_end_position"]
    info_dict["hit_in_idr"] = level_info_dict["hit_in_idr"]
    if level_info_dict["hit_in_idr"]:
        info_dict["idr_start"] = level_info_dict["idr_start"]
        info_dict["idr_end"] = level_info_dict["idr_end"]
    return info_dict


def ortholog_analysis_setup_with_odb_id(
    hits_df: pd.DataFrame,
    odb_main_folder,
    og_analysis_output_folder,
    root,
    logger,
    regex=None,
):
    odb_main_folder = Path(odb_main_folder)
    geneid_folders_dict = get_geneid_folder_dict(odb_main_folder)

    references_failed = {}
    references_passed = {}
    for i in hits_df.index:
        geneid = hits_df.loc[i, "odb_id"]
        uni_id = hits_df.loc[i, "UniprotID"]
        if not isinstance(uni_id, str):
            uni_id = "no_uniprot_id"
        ref_ind = int(hits_df.loc[i, "reference_index"])
        try:
            odb_folder = str(geneid_folders_dict[geneid])
        except KeyError as ke:
            logger.warning(f"no ortholog folder found for (gene_id: {geneid}) - {ke}")
            # hits_df.loc[i, 'odb_geneid'] = False
            references_failed[ref_ind] = "no ortholog folder found for gene id"
            continue

        hit_seq = hits_df.loc[i, "motif_match"]
        if "name" in hits_df.columns:
            name = hits_df.loc[i, "name"]
        else:
            name = geneid
        if "regex" in hits_df.columns:
            regex = hits_df.loc[i, "regex"]
        elif regex is not None:
            regex = regex
        ref_ind = int(hits_df.loc[i, "reference_index"])
        level_jsons = get_levels_json_dict(odb_folder)
        levels = get_sorted_levels_list(level_jsons)
        level_counts = get_number_seqs_in_all_levels(level_jsons)

        hit_output_folder = og_analysis_output_folder / f"{ref_ind}-{uni_id}_{geneid}"
        hit_output_folder.mkdir(exist_ok=True, parents=True)
        info_dict = {
            "UniprotID": uni_id,
            "query_gene_id": geneid,
            "odb_folder": odb_folder,
            "hit_sequence": hit_seq,
            "name": name,
            "reference_index": ref_ind,
            "regex": regex,
            "analysis_folder": str(hit_output_folder.resolve().relative_to(root)),
        }
        try:
            lvl_info_dict = create_level_info_dict(
                levels,
                level_jsons,
                hits_df.loc[i, "odb_mot_st"],
                hits_df.loc[i, "odb_mot_end"],
                **info_dict,
            )
        except IndexError as ie:
            logger.warning(f"index error for {geneid} - {ie}")
            references_failed[ref_ind] = str(ie)
            continue
        info_dict["level_info"] = lvl_info_dict
        info_dict["levels"] = levels
        info_dict["json_files_per_level"] = level_jsons
        info_dict["found hit"] = True
        info_dict = add_lvlinfo_2_info_dict(info_dict, lvl_info_dict[levels[0]])
        info_dict["num_sequences_per_level"] = level_counts
        group_info_json = hit_output_folder / "ortholog_group_info.json"
        with open(group_info_json, "w") as f:
            json.dump(info_dict, f, indent=4)

        references_passed[ref_ind] = str(group_info_json.resolve().relative_to(root))
    return references_failed, references_passed


def main(
    hits_file,
    og_analysis_output_folder,
    odb_og_run_folder,
    logger,
    meta_info: param_and_log_tools.Metainfo,
    root=here(),
    regex=None,
):
    meta_info.clear_meta_json()
    og_analysis_output_folder = Path(og_analysis_output_folder)
    og_analysis_output_folder.mkdir(exist_ok=True, parents=True)
    og_analysis_output_folder_rel_to_root = (
        og_analysis_output_folder.resolve().relative_to(root)
    )
    hits_df = import_and_reindex_hits_df(hits_file)
    (
        references_failed,
        references_passed,
    ) = ortholog_analysis_setup_with_odb_id(
        hits_df, odb_og_run_folder, og_analysis_output_folder, root, logger, regex=regex
    )
    meta_info.META_PARAMS["og_analysis_output_folder"] = str(
        og_analysis_output_folder_rel_to_root
    )
    meta_info.META_PARAMS["odb_main_folder"] = str(odb_og_run_folder)
    meta_info.META_PARAMS["ref_id_dict"] = references_passed
    meta_info.META_PARAMS["failed_id_dict"] = references_failed
    meta_info.META_PARAMS["notes"] = {}
    meta_info.save_meta_dict()


# if __name__ == "__main__":
# main(HITS_FILE, OG_ANALYSIS_OUTPUT_FOLDER, ODB_MAIN_FOLDER, logger)
# logger.debug("end of program")
