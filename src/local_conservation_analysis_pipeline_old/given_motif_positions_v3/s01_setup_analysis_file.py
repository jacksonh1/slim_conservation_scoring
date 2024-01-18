import json
from pathlib import Path

# import meta_tools
import pandas as pd
import param_and_log_tools
from pyprojroot import here

import local_env_variables.env_variables_and_filepaths as fp
import local_orthoDB_analysis_tools_v2.hit_sequence_tools as hit_tools


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
    info_dict["hit_in_idr"] = level_info_dict["hit_in_idr"]
    if level_info_dict["hit_in_idr"]:
        info_dict["idr_start"] = level_info_dict["idr_start"]
        info_dict["idr_end"] = level_info_dict["idr_end"]
    return info_dict


def create_analysis_dict(reference_index, odbid, hit_sequence, hit_start_position, hit_end_position, geneid_folders_dict, output_folder, uniprotid=None, name=None, regex=None):
    analysis_dict = {
        "reference_index": reference_index,
        "query_gene_id": odbid,
        "uniprotid": uniprotid,
        "name": name,
        "hit_sequence": hit_sequence,
        # "UniprotID": uniprotid,
        # "name": name,
        'hit_start_position': hit_start_position,
        'hit_end_position': hit_end_position,
        "failures": {},
    }
    output_filename = Path(output_folder) / 'analysis_jsons' / f"{reference_index}-{uniprotid}_{odbid.replace(':', '_')}.json"
    try:
        odb_folder = str(geneid_folders_dict[odbid])
    except KeyError as ke:
        analysis_dict["failures"]["critical"] = "no ortholog folder found for gene id"
        return analysis_dict, output_filename
    level_jsons = get_levels_json_dict(odb_folder)
    levels_list = get_sorted_levels_list(level_jsons)
    level_counts = get_number_seqs_in_all_levels(level_jsons)
    lvl_info_dict = create_level_info_dict(
        levels=levels_list,
        level_jsons_dict=level_jsons,
        hit_start_position=hit_start_position,
        hit_end_position=hit_end_position,
        **{k: v for k, v in analysis_dict.items() if k not in ['failures', 'hit_start_position', 'hit_end_position', 'hit_sequence']},
    )
    analysis_dict["level_info"] = lvl_info_dict
    analysis_dict["odb_group_info_files"] = level_jsons
    analysis_dict["num_sequences_per_level"] = level_counts
    analysis_dict["levels"] = levels_list
    analysis_dict = add_lvlinfo_2_info_dict(analysis_dict, lvl_info_dict)
    return analysis_dict, output_filename


def create_analysis_file(
    reference_index,
    odbid,
    hit_sequence,
    hit_start_position,
    hit_end_position,
    odb_og_run_folder,
    output_folder,
    root,
):
    geneid_folders_dict = get_geneid_folder_dict(odb_og_run_folder)
    output_folder = Path(output_folder)
    output_folder.mkdir(exist_ok=True, parents=True)
    analysis_dict, output_filename = create_analysis_dict(
        reference_index=reference_index,
        odbid=odbid,
        hit_sequence=hit_sequence,
        hit_start_position=hit_start_position,
        hit_end_position=hit_end_position,
        geneid_folders_dict=geneid_folders_dict,
        output_folder=output_folder
    )
    with open(output_filename, "w") as f:
        json.dump(analysis_dict, f, indent=4)
    output_file_relative_2_root = (
        output_filename.resolve().relative_to(root)
    )
    return output_file_relative_2_root




