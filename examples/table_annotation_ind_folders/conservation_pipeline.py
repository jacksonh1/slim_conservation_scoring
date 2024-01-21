import json
import multiprocessing
import shutil
from pathlib import Path

import local_config.pipeline_parameters as conf
import local_env_variables.env_variables as env
import s1setup_folder as setup
import s2define_idrs as find_idrs
import s3find_hit as find_hit
import s4add_lvlinfo as add_lvlinfo
import s5compute_scores as compute_scores
import s5multilevel_plots as multilevel_plots
import s6output_aln_slice as aln_slice
import s7map_results as map_results
import s8table_annotations as table_annotations
import yaml
from attrs import asdict

CONFIG_FILE = './params.yaml'
N_CORES = multiprocessing.cpu_count()


def load_config(config_file: str):
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


def run_pipeline(file, config: conf.PipelineParameters):
    # find idrs
    find_idrs.main(
        info_json_file=file,
        find_idrs=config.idr_params.find_idrs,
        idr_map_file=config.idr_params.idr_map_file,
        iupred_cutoff=config.idr_params.iupred_cutoff,
        gap_merge_threshold=config.idr_params.gap_merge_threshold,
        idr_min_length=config.idr_params.idr_min_length,
    )
    if config.hit_sequence_params.hit_sequence_search_method == "search":
        res = find_hit.main(
            json_file=file,
            longest_common_subsequence=config.hit_sequence_params.longest_common_subsequence,
            lcs_min_length=config.hit_sequence_params.lcs_min_length,
            target_hit_length=config.hit_sequence_params.target_hit_length,
        )
    else:
        res = "pass"
    if res == "fail":
        print(f"failed on {file}")
        return
    add_lvlinfo.main(file, config.filter_params.min_num_orthos)
    if len(config.score_methods) > 0:
        for scoremethod in config.score_methods:
            compute_scores.compute_conservation_scores(
                json_file=file,
                score_key=scoremethod.score_key,
                score_params=scoremethod.score_kwargs,
            )
    multilevel_plots.main(
        json_file=file,
        score_key=config.multilevel_plot_params.score_key,
    )
    aln_slice.main(
        json_file=file,
        n_flanking_cols=config.aln_slice_params.n_flanking_cols,
    )


def main(config_file, n_cores):
    config = load_config(config_file)
    if config.clear_files:
        if Path(config.output_folder).exists():
            shutil.rmtree(config.output_folder)

    setup.main(
        hits_file=config.table_file,
        database_key_file=config.database_filekey,
        output_folder=config.output_folder,
        hit_search_method=config.hit_sequence_params.hit_sequence_search_method,
    )
    json_files = get_passing_jsons(Path(config.output_folder))
    p = multiprocessing.Pool(n_cores)
    f_args = [(i, config) for i in json_files]
    p.starmap(run_pipeline, f_args)
    p.close()
    p.join()
    map_results.main(config.output_folder, config.table_file)
    table_annotations.main(
        config.table_file,
        config.table_annotation_params.score_key_for_table,
        config.table_annotation_params.levels,
    )


if __name__ == "__main__":
    main(CONFIG_FILE, N_CORES)