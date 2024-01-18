# %%
import json
import multiprocessing
import shutil
from pathlib import Path

from pyprojroot import here

import local_conservation_analysis_pipeline.given_motif_positions.param_and_log_tools as param_and_log_tools
import local_env_variables.env_variables_and_filepaths as fp
from local_orthoDB_analysis_tools_v2 import group_tools_v2 as group_tools


def copy_seq_file_2_analysis_folder(group_obj: group_tools.ortholog_group_info, level, sequence_file_key, output_folder):
    output_folder = Path(output_folder)
    level_info_o = group_obj.level_info_objects[level]
    src = level_info_o.sequence_files[sequence_file_key]
    dst = output_folder / f"{group_obj.query_gene_id.replace(':','_')}_{Path(src).name}"
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
        for seqkey in sequence_file_key:
            filename = copy_seq_file_2_analysis_folder(group_obj=og_obj, level=level, sequence_file_key=seqkey, output_folder=output_folder)
            group[seqkey] = str(filename.resolve().relative_to(root))
    og_obj._overwrite_json()


def main(
    levels_to_copy: str | list,
    sequence_file_key: list,
    og_analysis_output_folder: str|Path,
    meta_info: param_and_log_tools.Metainfo,
    root: Path,
    multiprocess=False,
):
    og_info_jsons = meta_info.jsons_to_analyze()
    og_info_jsons = [root / i for i in og_info_jsons]
    og_analysis_output_folder = Path(og_analysis_output_folder)
    alignment_folder = og_analysis_output_folder / "alignments"
    # og_analysis_output_folder_rel_to_root = (
    #     og_analysis_output_folder.resolve().relative_to(root)
    # )
    # alignment_folder = root / og_analysis_output_folder_rel_to_root / "alignments"
    alignment_folder.mkdir(exist_ok=True, parents=True)
    if multiprocess:
        p = multiprocessing.Pool(multiprocessing.cpu_count() - 1)
        f_args = [(i, alignment_folder, levels_to_copy, sequence_file_key, root) for i in og_info_jsons]
        p.starmap(copy_alignment_files, f_args)
        p.close()
        p.join()
    else:
        for og_info_json in og_info_jsons:
            print(f"Copying alignment files for {og_info_json}")
            copy_alignment_files(og_info_json, alignment_folder, levels_to_copy, sequence_file_key, root=root)
