import os
import sys
import warnings
from pathlib import Path

import local_conservation_analysis_pipeline_old.searching_for_hit_sequence.param_and_log_tools as param_and_log_tools
import local_conservation_analysis_pipeline_old.searching_for_hit_sequence.s01_v2_setup_analysis as setup_analysis
import local_conservation_analysis_pipeline_old.searching_for_hit_sequence.s01_v4_setup_lvl_hit_info as setup_hit_info
import local_conservation_analysis_pipeline_old.searching_for_hit_sequence.s02_filter as s02filter
import local_conservation_analysis_pipeline_old.searching_for_hit_sequence.s03_score_fl_alns as score_fl_alns_valday
import local_conservation_analysis_pipeline_old.searching_for_hit_sequence.s03_score_fl_alns_PE as score_fl_alns_PE
import local_conservation_analysis_pipeline_old.searching_for_hit_sequence.s04_make_multi_level_plots as multi_level_plots
import yaml
import local_env_variables.env_variables as env

# sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


# warnings.filterwarnings("ignore", module = "matplotlib\..*" )


ROOT = env.ROOT
params = param_and_log_tools.load_params("./params.yaml")
PARAM_OBJ = param_and_log_tools.RunParameters_dict(params)
META_INFO = param_and_log_tools.Metainfo(PARAM_OBJ.META_FILE)


def main(meta_info, param_obj, root):
    if "s01" in param_obj.PARAMS["scripts_to_run"]:
        output_folder = Path(param_obj.OG_ANALYSIS_OUTPUT_FOLDER) / "orthogroup_analysis_jsons"
        setup_analysis.main(
            hits_file=param_obj.HITS_FILE,
            orthogroup_analysis_json_folder=output_folder,
            odb_og_run_folder=param_obj.ODB_MAIN_FOLDER,
            logger=logger,
            meta_info=meta_info,
            root=root,
        )
        meta_info = meta_info.reload_meta_dict()
    if "s02" in param_obj.PARAMS["scripts_to_run"]:
        s02filter.main(
            meta_info,
            root=root,
        )
        meta_info = meta_info.reload_meta_dict()
    if "s03" in param_obj.PARAMS["scripts_to_run"]:
        score_fl_alns_valday.main(
            matrix_name=param_obj.SCORE_PARAMS["matrix_name"],
            gap_frac_cutoff=param_obj.SCORE_PARAMS["gap_frac_cutoff"],
            output_folder=Path(param_obj.OG_ANALYSIS_OUTPUT_FOLDER) / f"conservation_scores-{param_obj.SCORE_PARAMS['matrix_name']}-g_{param_obj.SCORE_PARAMS['gap_frac_cutoff']}",
            root=root,
            meta_info=meta_info,
            multiprocess=True,
        )
        score_fl_alns_PE.main(
            gap_frac_cutoff=param_obj.SCORE_PARAMS["gap_frac_cutoff"],
            output_folder=Path(param_obj.OG_ANALYSIS_OUTPUT_FOLDER) / f"conservation_scores-property_entropy-g_{param_obj.SCORE_PARAMS['gap_frac_cutoff']}",
            root=root,
            meta_info=meta_info,
            multiprocess=True,
        )
    if "s04" in param_obj.PARAMS["scripts_to_run"]:
        multi_level_plots.main(meta_info=meta_info, root=root)
    if "s05" in param_obj.PARAMS["scripts_to_run"]:
        add_best_level.main(
            param_obj.REQUIRED_USABLE_FRACTION, meta_info=meta_info, root=root
        )
        meta_info = meta_info.reload_meta_dict()
    if "s06" in param_obj.PARAMS["scripts_to_run"]:
        copy_alignments.main(
            levels_to_copy=param_obj.LEVELS_TO_COPY,
            sequence_file_key=param_obj.SEQUENCE_FILE_KEY,
            alignment_folder=Path(param_obj.OG_ANALYSIS_OUTPUT_FOLDER) / 'alignments',
            root=root,
            meta_info=meta_info,
            multiprocess=False,
        )


if __name__ == "__main__":
    main(META_INFO, PARAM_OBJ, ROOT)


# @hydra.main(config_path="config.yaml")
# def pipeline(cfg):
#     print(cfg.pretty())
