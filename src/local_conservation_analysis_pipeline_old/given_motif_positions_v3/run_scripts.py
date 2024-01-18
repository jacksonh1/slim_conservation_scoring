import warnings

import multiprocessing
import s01_v2_setup_analysis as setup_analysis
import s02_copy_over_alignments as copy_alignments
import s03_filter as s03filter
import s04_score_fl_alns as score_fl_alns_valday
import s04_score_fl_alns_PE as score_fl_alns_PE
import s05_make_multi_level_plots as multi_level_plots
import yaml

warnings.filterwarnings("ignore", module = "matplotlib\..*" )

import pipeline_tools
from local_env_variables import env_variables_and_filepaths as fp

ROOT = fp.conservation_analysis_root
params = pipeline_tools.load_params('./params.yaml')
param_obj = pipeline_tools.RunParameters_dict(params)
meta_info = pipeline_tools.Metainfo(param_obj.META_FILE)

def remove_failed_hits_from_meta_dict(
    og_info_json,
    meta_info: pipeline_tools.Metainfo
):
    with open(og_info_json, 'r') as f:
        og_info = json.load(f)
    if len(og_info['levels']) == 0:
        meta_info.add_id_to_fail_dict(og_info["reference_index"], f'no levels with >=10 orthologs')
    if not og_info['hit_in_idr']:
        meta_info.add_id_to_fail_dict(og_info["reference_index"], 'hit sequence is not in an IDR')
    return meta_info

def run_01(param_obj, meta_info, root):
    if 's01' in param_obj.PARAMS['scripts_to_run']:
        logger = pipeline_tools.setup_logging("./logs/s01.log")
        meta_info = setup_analysis.main(
            hits_file=param_obj.HITS_FILE,
            og_analysis_output_folder=param_obj.OG_ANALYSIS_OUTPUT_FOLDER,
            odb_og_run_folder=param_obj.ODB_MAIN_FOLDER,
            logger=logger,
            meta_info=meta_info,
            root=root,
            regex=param_obj.PARAMS["regex"],
        )
    return meta_info

def non_multiprocessed(og_info_json, param_obj, meta_info, root):
    if 's02' in param_obj.PARAMS['scripts_to_run']:
        logger = pipeline_tools.setup_logging("./logs/s02.log")
        copy_alignments.main(
            levels_to_copy=param_obj.LEVELS_TO_COPY,
            sequence_file_key=param_obj.SEQUENCE_FILE_KEY,
            og_info_json=og_info_json,
            root=root,
        )
    if 's03' in param_obj.PARAMS['scripts_to_run']:
        s03filter.main(
            og_info_json,
            root=root,
        )
    meta_info = remove_failed_hits_from_meta_dict(og_info_json, meta_info)
    return meta_info


def multiprocessed(og_info_json, param_obj, root):
    if 's04' in param_obj.PARAMS['scripts_to_run']:
        score_fl_alns_valday.main(
            matrix_name=param_obj.SCORE_PARAMS["matrix_name"],
            gap_frac_cutoff=param_obj.SCORE_PARAMS["gap_frac_cutoff"],
            og_info_json=og_info_json,
            root=root,
        )
        score_fl_alns_PE.main(
            gap_frac_cutoff=param_obj.SCORE_PARAMS["gap_frac_cutoff"],
            og_info_json=og_info_json,
            root=root,
        )
    if 's05' in param_obj.PARAMS['scripts_to_run']:
        multi_level_plots.main(og_info_json, root=root)


def main(param_obj, meta_info, root):
    meta_info = run_01(param_obj, meta_info, root)
    og_info_jsons = [root / i for i in meta_info.jsons_to_analyze()]
    for og_info_json in og_info_jsons:
        meta_info = non_multiprocessed(og_info_json, param_obj, meta_info, root)
    og_info_jsons = [root / i for i in meta_info.jsons_to_analyze()]
    with multiprocessing.Pool(16) as p:
        f_args = [(i, param_obj, root) for i in og_info_jsons]
        results_async = p.starmap_async(multiprocessed, f_args)
        for result in results_async.get():
            pass

if __name__ == "__main__":
    main(param_obj, meta_info, ROOT)

# @hydra.main(config_path="config.yaml")
# def pipeline(cfg):
#     print(cfg.pretty())
