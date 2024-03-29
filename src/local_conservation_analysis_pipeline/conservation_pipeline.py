import argparse
import json
import multiprocessing
import shutil
from functools import partial
from pathlib import Path

import yaml
from attrs import asdict

import local_config.conservation_pipeline_parameters as conf
from local_conservation_analysis_pipeline import (s1setup_folder,
                                                  s2define_idrs, s3find_hit,
                                                  s4add_lvlinfo,
                                                  s5compute_scores,
                                                  s6multilevel_plots,
                                                  s7output_aln_slice,
                                                  s8calculate_annotations,
                                                  s9add_annotations2table)

CONFIG_FILE = "./params.yaml"
N_CORES = multiprocessing.cpu_count()


def load_config(config_file: str) -> conf.PipelineParameters:
    # if config_file is None:
    #     config = conf.PipelineParameters()
    # else:
    # with open(config_file, 'r') as f:
    #     config_dict = yaml.safe_load(f)
    # config = conf.PipelineParameters.from_dict(config_dict)
    with open(config_file, "r") as f:
        config_dict = yaml.safe_load(f)
    config = conf.PipelineParameters.from_dict(config_dict)
    return config


def get_passing_jsons(search_dir):
    search_dir = Path(search_dir)
    json_files = search_dir.rglob("*.json")
    passing_jsons = []
    for json_file in json_files:
        with open(json_file, "r") as f:
            json_dict = json.load(f)
        if "critical_error" not in json_dict:
            if "reference_index" in json_dict:
                passing_jsons.append(json_file)
    return passing_jsons


def run_multiprocess_steps(file, config: conf.PipelineParameters):
    # find idrs
    if "s2" in config.steps_to_run:
        s2define_idrs.main(
            info_json_file=file,
            find_idrs=config.idr_params.find_idrs,
            idr_map_file=config.idr_params.idr_map_file,
            iupred_cutoff=config.idr_params.iupred_cutoff,
            gap_merge_threshold=config.idr_params.gap_merge_threshold,
            idr_min_length=config.idr_params.idr_min_length,
        )
    if "s3" in config.steps_to_run:
        res = s3find_hit.main(
            json_file=file,
            search_method=config.hit_sequence_params.hit_sequence_search_method,
            longest_common_subsequence=config.hit_sequence_params.longest_common_subsequence,
            lcs_min_length=config.hit_sequence_params.lcs_min_length,
            target_hit_length=config.hit_sequence_params.target_hit_length,
        )
        if res == "fail":
            print(f"failed on {file}")
            return
    if "s4" in config.steps_to_run:
        s4add_lvlinfo.main(file, config.filter_params.min_num_orthos)
    if "s5" in config.steps_to_run:
        if len(config.score_methods) > 0:
            for scoremethod in config.score_methods:
                s5compute_scores.compute_conservation_scores(
                    json_file=file,
                    score_key=scoremethod.score_key,
                    score_params=scoremethod.score_kwargs,
                )
    if "s6" in config.steps_to_run:
        s6multilevel_plots.multi_level_plots(
            json_file=file,
            score_key=config.multilevel_plot_params.score_key,
            score_type=config.multilevel_plot_params.score_type,
            num_bg_scores_cutoff=config.multilevel_plot_params.num_bg_scores_cutoff,
        )
    if "s7" in config.steps_to_run:
        s7output_aln_slice.main(
            json_file=file,
            n_flanking_aas=config.aln_slice_params.n_flanking_aas,
            whole_idr=config.aln_slice_params.whole_idr,
        )


def main(config_file, n_cores):
    config = load_config(config_file)

    if "s1" in config.steps_to_run:
        if config.clear_files:
            if Path(config.output_folder).exists():
                shutil.rmtree(config.output_folder)
        s1setup_folder.main(
            hits_file=config.table_file,
            database_key_file=config.database_filekey,
            output_folder=config.output_folder,
            hit_search_method=config.hit_sequence_params.hit_sequence_search_method,
        )
    json_files = get_passing_jsons(Path(config.output_folder))
    # p = multiprocessing.Pool(n_cores)
    # f_args = [(i, config) for i in json_files]
    # p.starmap(run_multiprocess_steps, f_args)
    # p.close()
    # p.join()
    with multiprocessing.Pool(n_cores) as p:
        results_iterator = p.imap_unordered(
            partial(run_multiprocess_steps, config=config), json_files, chunksize=1
        )
        for result in results_iterator:
            pass
    if "s8" in config.steps_to_run:
        print("calculating annotations")
        s8calculate_annotations.main(
            main_output_folder=config.output_folder,
            image_score_key=config.multilevel_plot_params.score_key,
            table_annotation_score_key=config.table_annotation_params.score_key_for_table,
            regex=config.table_annotation_params.motif_regex,
        )
    if "s9" in config.steps_to_run:
        s9add_annotations2table.main(
            annotations_file=Path(config.output_folder) / "annotations.json",
            table_file=config.table_file,
            table_annotations=config.table_annotation_params.annotations,
            table_annotation_levels=config.table_annotation_params.levels,
            output_table_file=str(config.table_file).replace(".csv", "_ANNOTATED.csv"),
        )
    if config.clean_analysis_files:
        shutil.rmtree(config.output_folder)
    # delete the reindexed table file
    Path(str(config.table_file).replace(".csv", "_original_reindexed.csv")).unlink()
    # save config to a parameters file
    with open(Path(config.output_folder) / 'processing_parameters.yaml', 'w') as f:
        yaml.dump(asdict(config), f, default_flow_style=False)


if __name__ == "__main__":
    main(CONFIG_FILE, N_CORES)
