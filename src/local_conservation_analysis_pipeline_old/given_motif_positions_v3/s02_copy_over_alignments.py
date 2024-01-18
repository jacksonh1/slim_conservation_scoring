# %%
import json
import multiprocessing
import shutil
from pathlib import Path

import pipeline_tools
from pyprojroot import here

import local_env_variables.env_variables_and_filepaths as fp
from local_orthoDB_analysis_tools_v2 import group_tools_v2 as group_tools

# def copy_seq_file_2_analysis_folder(src_file, dst_folder):
#     src = src_file
#     dst = Path(dst_folder) / Path(src).name
#     shutil.copy(src, dst)
#     return dst

def copy_seq_file_2_analysis_folder(level_info_o: group_tools.level_info, sequence_file_key, aln_folder):
    src = level_info_o.sequence_files[sequence_file_key]
    dst = Path(aln_folder) / f"{level_info_o.reference_index}-{level_info_o.query_gene_id.replace(':','_')}_{Path(src).name}"
    shutil.copy(src, dst)
    return dst

def copy_alignment_files(
    og_info_json,
    output_folder,
    levels: str | list = "all",
    sequence_file_key=["fasta alignment - OG LDO cdhit"],
    root=here(),
):
    og_obj = group_tools.ortholog_group_info(og_info_json, root=root)
    og_obj.load_all_level_info_objects()
    if isinstance(levels, str) and levels != "all":
        if levels in og_obj.levels:
            levels = [levels]
        else:
            raise ValueError(f"{levels} is not a valid level")
    # assert levels == 'all' or isinstance(levels, list)
    if levels == "all":
        levels = og_obj.levels
    elif isinstance(levels, list):
        assert all(
            [level in og_obj.levels for level in levels]
        ), f"Not all levels are valid {levels}"
    if isinstance(sequence_file_key, str):
        sequence_file_key = [sequence_file_key]

    output_folder = Path(output_folder)
    aln_folder = output_folder / "alignments"
    for level in levels:
        group = og_obj.info_dict["level_info"][level].setdefault(
            "local_sequence_files", {}
        )
        level_info_o = og_obj.level_info_objects[level]
        for seqkey in sequence_file_key:
            filename = copy_seq_file_2_analysis_folder(level_info_o, seqkey, aln_folder)
            group[seqkey] = str(filename.resolve().relative_to(root))
    og_obj._overwrite_json()


def main(
    levels_to_copy: str | list,
    sequence_file_key: list,
    og_info_json: str|Path,
    root: Path,
):
    og_info_json = root / og_info_json
    copy_alignment_files(og_info_json, levels_to_copy, sequence_file_key, root)

