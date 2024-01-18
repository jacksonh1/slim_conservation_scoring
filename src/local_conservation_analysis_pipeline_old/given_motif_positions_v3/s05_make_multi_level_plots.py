import json
import multiprocessing
import os
import re
import subprocess
import sys
from pathlib import Path

import local_env_variables.env_variables_and_filepaths as fp
import local_seqtools.general_utils as tools
import matplotlib.pyplot as plt
import pipeline_tools
from Bio import Align, AlignIO, Seq, SeqIO
from local_orthoDB_analysis_tools_v2 import \
    conservation_score_plots as cons_plots
from local_orthoDB_analysis_tools_v2 import group_tools_v2 as group_tools
from pyprojroot import here

ROOT = fp.conservation_analysis_root


def multi_level_driver(ortholog_group_info_json, root):
    SCORE_KEY='property_entropy_score_json'
    og_obj = group_tools.ortholog_group_info(ortholog_group_info_json, root=ROOT)
    if not og_obj.found_hit:
        print(f'No hit found for {ortholog_group_info_json}')
        return
    og_obj.load_all_level_score_objects_with_z_score(score_key=SCORE_KEY)

    multi_plot_filename_1 = cons_plots.ogscreen_plot_from_og_folder_class(og_obj, score_name=f'{SCORE_KEY.replace("_score_json","")}-stripped', bar_ylim=[-0.5, 2.5], strip_gaps=True)
    og_obj.add_item_to_json('multi_plot_filename_1', str(multi_plot_filename_1.resolve().relative_to(ROOT)), save_json=True)

    multi_plot_filename_2 = cons_plots.ogscreen_plot_from_og_folder_class(og_obj, score_name=f'{SCORE_KEY.replace("_score_json","")}', bar_ylim=[-0.5, 2.5], strip_gaps=False)
    og_obj.add_item_to_json('multi_plot_filename_2', str(multi_plot_filename_2.resolve().relative_to(ROOT)), save_json=True)

    SCORE_KEY='asym_valday_score_json'
    og_obj = group_tools.ortholog_group_info(ortholog_group_info_json, root=ROOT)
    if not og_obj.found_hit:
        print(f'No hit found for {ortholog_group_info_json}')
        return
    og_obj.load_all_level_score_objects_with_z_score(score_key=SCORE_KEY)

    
    multi_plot_filename_3 = cons_plots.ogscreen_plot_from_og_folder_class(og_obj, score_name=f'{SCORE_KEY.replace("_score_json","")}-stripped', bar_ylim=[-0.5, 2.5], strip_gaps=True)
    og_obj.add_item_to_json('multi_plot_filename_3', str(multi_plot_filename_3.resolve().relative_to(ROOT)), save_json=True)

    multi_plot_filename_4 = cons_plots.ogscreen_plot_from_og_folder_class(og_obj, score_name=f'{SCORE_KEY.replace("_score_json","")}', bar_ylim=[-0.5, 2.5], strip_gaps=False)
    og_obj.add_item_to_json('multi_plot_filename_4', str(multi_plot_filename_4.resolve().relative_to(ROOT)), save_json=True)
    plt.close("all")



def main(
    og_info_json: str|Path,
    root: Path,
):
    og_info_json = root / og_info_json
    multi_level_driver(og_info_json, root=root)
