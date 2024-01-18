# %%
import multiprocessing
import shutil
from pathlib import Path

from pyprojroot import here

import local_conservation_analysis_pipeline.searching_for_hit_sequence.param_and_log_tools as param_and_log_tools
from local_orthoDB_analysis_tools_v2 import group_tools_v2 as group_tools

# def copy_seq_file_2_analysis_folder(src_file, dst_folder):
#     src = src_file
#     dst = Path(dst_folder) / Path(src).name
#     shutil.copy(src, dst)
#     return dst

def copy_seq_file_2_analysis_folder(level_info_o, sequence_file_key):
    src = level_info_o.sequence_files[sequence_file_key]
    dst = level_info_o.analysis_folder / Path(src).name
    shutil.copy(src, dst)
    return dst

def copy_alignment_files(
    og_info_json,
    levels: str | list = "all",
    sequence_file_key=["fasta alignment - OG LDO cdhit"],
    root=here(),
):
    og_obj = group_tools.ortholog_group_info(og_info_json, root=root)
    og_obj.load_all_level_info_objects()
    if isinstance(levels, str) and levels != "all" and levels != 'best':
        if levels in og_obj.levels:
            levels = [levels]
        else:
            raise ValueError(f"{levels} is not a valid level")
    # assert levels == 'all' or isinstance(levels, list)
    if levels == "all":
        levels = og_obj.levels
    elif levels == 'best':
        if og_obj.best_level is None:
            print(f"No best level for {og_obj}")
            return
        levels = [og_obj.best_level]
    elif isinstance(levels, list):
        assert all(
            [level in og_obj.levels for level in levels]
        ), f"Not all levels are valid {levels}"
    if isinstance(sequence_file_key, str):
        sequence_file_key = [sequence_file_key]

    for level in levels:
        group = og_obj.info_dict["level_info"][level].setdefault(
            "local_sequence_files", {}
        )
        level_info_o = og_obj.level_info_objects[level]
        for seqkey in sequence_file_key:
            filename = copy_seq_file_2_analysis_folder(level_info_o, seqkey)
            group[seqkey] = str(filename.resolve().relative_to(root))
    og_obj._overwrite_json()


def main(
    levels_to_copy: str | list,
    sequence_file_key: list,
    meta_info: param_and_log_tools.Metainfo,
    root: Path,
    multiprocess=False,
):
    og_info_jsons = meta_info.jsons_to_analyze()
    og_info_jsons = [root / i for i in og_info_jsons]
    if multiprocess:
        p = multiprocessing.Pool(multiprocessing.cpu_count() - 1)
        f_args = [(i, levels_to_copy, sequence_file_key, root) for i in og_info_jsons]
        p.starmap(copy_alignment_files, f_args)
        p.close()
        p.join()
    else:
        for og_info_json in og_info_jsons:
            print(f"Copying alignment files for {og_info_json}")
            copy_alignment_files(og_info_json, levels_to_copy, sequence_file_key, root=root)
