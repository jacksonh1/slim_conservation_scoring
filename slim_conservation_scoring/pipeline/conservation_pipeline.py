import json
import multiprocessing
import shutil
from functools import partial
from pathlib import Path

import yaml
from attrs import asdict
import copy
from typing import Callable

import slim_conservation_scoring.config.conservation_pipeline_parameters as conf
from slim_conservation_scoring.pipeline import (
    s1setup_folder,
    s2define_idrs,
    s3find_hit,
    s4add_lvlinfo,
    s5a_compute_aln_scores,
    s5b_run_pairk_aln,
    s5c_add_hit_scores,
    s6multilevel_plots,
    s7output_aln_slice,
    s8calculate_annotations,
    s9add_annotations2table,
)
import time

# import conservation_scores.tools.esm_model as esm_model
from slim_conservation_scoring.seqtools import esm_tools

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


def get_passing_jsons(search_dir, exclude_dir: str | Path | None = None):
    search_dir = Path(search_dir)
    json_files = search_dir.rglob("*.json")
    passing_jsons = []
    for json_file in json_files:
        if "clustered" in json_file.name:
            continue
        if exclude_dir is not None:
            if str(exclude_dir) in str(json_file):
                continue
        with open(json_file, "r") as f:
            json_dict = json.load(f)
        if "critical_error" not in json_dict:
            if "reference_index" in json_dict:
                passing_jsons.append(json_file)
    return passing_jsons


def remove_failed_jsons(json_files):
    passing_jsons = []
    for json_file in json_files:
        with open(json_file, "r") as f:
            json_dict = json.load(f)
        if "critical_error" not in json_dict:
            if (
                "reference_index" in json_dict
                and len(json_dict["levels_passing_filters"]) > 0
            ):
                passing_jsons.append(json_file)
    return passing_jsons


def run_setup_steps(file, config: conf.PipelineParameters):
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


def calculate_aln_cons_scores(file, config: conf.PipelineParameters, **kwargs):
    if len(config.msa_score_methods) > 0:
        for scoremethod in config.msa_score_methods:
            s5a_compute_aln_scores.compute_aln_cons_scores(
                json_file=file,
                score_key=scoremethod.score_key,
                function_name=scoremethod.function_name,
                function_params=scoremethod.function_params,
                level=scoremethod.level,
                **kwargs,
            )


def run_pairk_aln(file, config: conf.PipelineParameters, **kwargs):
    if len(config.pairk_aln_methods) > 0:
        for pairk_method in config.pairk_aln_methods:
            s5b_run_pairk_aln.pairwise_kmer_alignment_driver(
                json_file=file,
                score_key=pairk_method.score_key,
                function_name=pairk_method.function_name,
                function_params=pairk_method.function_params,
                level=pairk_method.level,
                lflank=pairk_method.lflank,
                rflank=pairk_method.rflank,
                **kwargs,
            )


def run_pairk_embedding_aln(file, config: conf.PipelineParameters, **kwargs):
    if "embedding_file_path" in kwargs:
        if kwargs["embedding_file_path"] is not None:
            mod = None
        else:
            mod = esm_tools.ESM_Model(
                model_name=config.esm_params.model_name,
                threads=config.esm_params.threads,
            )
    else:
        mod = esm_tools.ESM_Model(
            model_name=config.esm_params.model_name, threads=config.esm_params.threads
        )
    score_kwargs = dict(
        device=config.esm_params.device,
        mod=mod,
    )
    for scoremethod in config.pairk_embedding_aln_methods:
        score_kwargs.update(scoremethod.function_params)
        s5b_run_pairk_aln.pairwise_kmer_alignment_driver(
            json_file=file,
            score_key=scoremethod.score_key,
            function_name=scoremethod.function_name,
            function_params=score_kwargs,
            level=scoremethod.level,
            lflank=scoremethod.lflank,
            rflank=scoremethod.rflank,
            **kwargs,
        )


def calculate_pairk_conservation(file, config: conf.PipelineParameters):
    if len(config.pairk_aln_methods) > 0:
        for pairk_method in config.pairk_aln_methods:
            s5c_add_hit_scores.compute_hit_conservation_scores(
                json_file=file,
                pairk_method=pairk_method,
                params=config.pairk_conservation_params,
            )
    if len(config.pairk_embedding_aln_methods) > 0:
        for pairk_method in config.pairk_embedding_aln_methods:
            s5c_add_hit_scores.compute_hit_conservation_scores(
                json_file=file,
                pairk_method=pairk_method,
                params=config.pairk_conservation_params,
            )


def generate_plots(file, config: conf.PipelineParameters):
    if "s6" in config.steps_to_run:
        if config.multilevel_plot_params.score_key is not None:
            s6multilevel_plots.multi_level_plot_driver(
                json_file=file,
                config=config,
            )
        else:
            print("no score key specified for multilevel plots")
    if "s7" in config.steps_to_run:
        s7output_aln_slice.main(
            json_file=file,
            n_flanking_aas=config.aln_slice_params.n_flanking_aas,
            whole_idr=config.aln_slice_params.whole_idr,
        )


def multiprocess_function(
    func: Callable,
    config: conf.PipelineParameters,
    json_files: list,
    n_cores: int,
    *args,
    **kwargs,
):
    with multiprocessing.Pool(n_cores) as p:
        results_iterator = p.imap_unordered(
            partial(func, config=config, *args, **kwargs), json_files, chunksize=1
        )
        for result in results_iterator:
            pass


def step1(config: conf.PipelineParameters, reindexed_table_file: Path):
    if config.new_index and not config.clear_files:
        raise ValueError(
            "new_index is True, but clear_files is False. If you want to reindex the table, you must set clear_files to True. Otherwise, you could have multiple reference indexes and all sorts of problems."
        )
    if config.clear_files:
        if Path(config.output_folder).exists():
            shutil.rmtree(config.output_folder)
            # shutil.rmtree(score_output_folder)
    if config.new_index:
        print("setting up folders and reindexing table file")
        s1setup_folder.main(
            hits_file=config.table_file,
            database_key_file=config.database_filekey,
            output_folder=config.output_folder,
            hit_search_method=config.hit_sequence_params.hit_sequence_search_method,
            new_index=config.new_index,
        )
    elif not config.new_index and not reindexed_table_file.exists():
        # this is for if your table already has a reference_index column that you want to use
        print(
            f"new_index is {config.new_index}, and a reindexed table file does not exist: {reindexed_table_file}, so the input table file should already have a reference_index column."
        )
        print("using the input table")
        s1setup_folder.main(
            hits_file=config.table_file,
            database_key_file=config.database_filekey,
            output_folder=config.output_folder,
            hit_search_method=config.hit_sequence_params.hit_sequence_search_method,
            new_index=config.new_index,
        )
        # reindexed_table_file = Path(config.table_file)
    elif not config.new_index and reindexed_table_file.exists():
        print(
            f"new_index is {config.new_index}, and the reindexed table file already exists: {reindexed_table_file}"
        )
        print("using the existing reindexed table file")
        s1setup_folder.main(
            hits_file=reindexed_table_file,
            database_key_file=config.database_filekey,
            output_folder=config.output_folder,
            hit_search_method=config.hit_sequence_params.hit_sequence_search_method,
            new_index=config.new_index,
        )
    # return reindexed_table_file


def main(config_file, n_cores):
    a = time.time()
    config = load_config(config_file)
    config.print_params()
    table_filename = Path(config.table_file).name
    reindexed_table_file = Path(config.output_folder) / table_filename.replace(
        ".csv", "_original_reindexed.csv"
    )
    annotated_table_file = Path(config.output_folder).parent / table_filename.replace(
        ".csv", "_ANNOTATED.csv"
    )
    Path(config.output_folder).mkdir(exist_ok=True, parents=True)
    score_output_folder = Path(config.output_folder) / "conservation_score_files"
    score_output_folder.mkdir(exist_ok=True)
    if "s1" in config.steps_to_run:
        step1(config, reindexed_table_file)
    json_files = get_passing_jsons(
        Path(config.output_folder), exclude_dir=score_output_folder
    )
    multiprocess_function(run_setup_steps, config, json_files, n_cores)
    # need to remove some more jsons that failed in `run_setup_steps`
    json_files = remove_failed_jsons(json_files)
    if "s5" in config.steps_to_run:
        print("calculating MSA scores")
        multiprocess_function(
            calculate_aln_cons_scores,
            config,
            json_files,
            n_cores,
            score_output_folder=score_output_folder,
        )
        print("running pairk alignment")
        multiprocess_function(
            run_pairk_aln,
            config,
            json_files,
            n_cores,
            score_output_folder=score_output_folder,
        )

    if len(config.pairk_embedding_aln_methods) > 0 and "s5" in config.steps_to_run:
        print("running pairk alignment with residue embeddings")
        multiprocess_function(
            run_pairk_embedding_aln,
            config,
            json_files,
            n_cores=config.esm_params.processes,
            score_output_folder=score_output_folder,
        )

    if "s5" in config.steps_to_run:
        print("calculating kmer scores from pairk alignments")
        multiprocess_function(calculate_pairk_conservation, config, json_files, n_cores)

    json_files = remove_failed_jsons(json_files)
    multiprocess_function(generate_plots, config, json_files, n_cores)

    if "s8" in config.steps_to_run:
        print("calculating annotations")
        s8calculate_annotations.main(
            main_output_folder=config.output_folder,
            image_score_key=config.multilevel_plot_params.score_key,
            table_annotation_score_key=config.table_annotation_params.score_key_for_table,
            regex=config.table_annotation_params.motif_regex,
        )
    if "s9" in config.steps_to_run:
        print("adding annotations to table")
        if not reindexed_table_file.exists():
            input_annotation_table = Path(config.table_file)
        else:
            input_annotation_table = reindexed_table_file
        s9add_annotations2table.main(
            annotations_file=Path(config.output_folder) / "annotations.json",
            table_file=input_annotation_table,
            table_annotations=config.table_annotation_params.annotations,
            table_annotation_levels=config.table_annotation_params.levels,
            output_table_file=annotated_table_file,
        )
    if config.clean_analysis_files:
        shutil.rmtree(config.output_folder)
    # delete the reindexed table file
    # if reindexed_table_file.exists():
    #     reindexed_table_file.unlink()
    # save config to a parameters file
    shutil.copyfile(
        config_file, Path(config.output_folder) / "processing_parameters_user.yaml"
    )
    with open(Path(config.output_folder) / "processing_parameters_full.yaml", "w") as f:
        x = asdict(config)
        yaml.dump(x, f)

    print("done")
    b = time.time()
    time_elapsed = b - a
    print(f"total time elapsed: {time_elapsed/60} minutes ({time_elapsed} seconds)")


if __name__ == "__main__":
    main(CONFIG_FILE, N_CORES)
