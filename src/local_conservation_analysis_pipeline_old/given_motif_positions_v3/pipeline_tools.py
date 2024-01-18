'''
================IN DEVELOPMENT==================
handle args in scripts

'''

import copy
import json
import logging
import os
from pathlib import Path

import yaml

DEFAULT_PARAMS = {
    "meta_file":"ortholog_analysis_meta.json",
    "hits_file":"",
    "odb_main_folder":"/home/jch/Documents/SLiM_bioinformatics_workspace/data/orthoDB_ortholog_group_runs/2023-07-19-whole_proteome/orthoDB_analysis_multiprocessed_minfrac-0.75/msa_by_organism",
    "og_analysis_output_folder":"ortholog_analysis",
    "hit_params":{
        "target_hit_length": 0
    },
    "score_params":{
        "matrix_name": "EDSSMat50_max_off_diagonal_norm",
        "gap_frac_cutoff":0.2
    },
    "required_usable_fraction":0.5,
    "sequence_file_key": [
        "fasta alignment - OG LDO cdhit"
    ],
    "levels_to_copy":"Vertebrata",
    "level_to_use_for_scores":"Vertebrata",
    "score_key_for_table":"asym_valday_score_json"
}


def load_params(user_parameters, default_parameters=DEFAULT_PARAMS):
    if isinstance(user_parameters, str):
        with open(user_parameters, "r") as f:
            user_parameters = yaml.safe_load(f)
    
    for k,v in default_parameters.items():
        if k not in user_parameters.keys():
            print(f"`{k}` parameter not provided")
            print(f"\t- using default value: `{v}`\n")

    # update default params with params provided by the user
    params = default_parameters.copy()
    params.update(user_parameters)
    return params


def setup_logging(log_file, level=logging.WARNING):
    log_file_path = Path(log_file)
    log_file_path.parent.mkdir(exist_ok=True, parents=True)
    logging.basicConfig(filename=log_file_path, level=level, 
                        format='%(asctime)s %(levelname)s %(name)s %(message)s')
    # clear the log file
    open(log_file_path, 'w').close()
    logging.debug('Start of program')
    logger=logging.getLogger(__name__)
    return logger


class RunParameters_dict:

    def __init__(self, parameter_dict):
        self.PARAMS = parameter_dict
        self.HITS_FILE = parameter_dict['hits_file']
        self.ODB_MAIN_FOLDER = parameter_dict['odb_main_folder']
        self.OG_ANALYSIS_OUTPUT_FOLDER = parameter_dict['og_analysis_output_folder']
        self.META_FILE = parameter_dict['meta_file']
        self.SCORE_PARAMS = parameter_dict['score_params']
        self.HIT_PARAMS = parameter_dict['hit_params']
        self.REQUIRED_USABLE_FRACTION = float(parameter_dict['required_usable_fraction'])
        self.SEQUENCE_FILE_KEY = parameter_dict['sequence_file_key']
        if isinstance(self.SEQUENCE_FILE_KEY, str):
            self.SEQUENCE_FILE_KEY = [self.SEQUENCE_FILE_KEY]
        self.LEVELS_TO_COPY = parameter_dict['levels_to_copy']
        self.LEVEL_TO_USE_FOR_SCORES = parameter_dict['level_to_use_for_scores']
        self.SCORE_KEY_FOR_TABLE = parameter_dict['score_key_for_table']



class RunParameters_yaml:

    def __init__(self, parameter_yaml_file):
        self.PARAMETER_yaml_FILE = parameter_yaml_file
        with open(parameter_yaml_file) as f:
            params = yaml.safe_load(f)
        self.PARAMS = params
        self.HITS_FILE = params['hits_file']
        self.ODB_MAIN_FOLDER = params['odb_main_folder']
        self.OG_ANALYSIS_OUTPUT_FOLDER = params['og_analysis_output_folder']
        self.META_FILE = params['meta_file']
        self.SCORE_PARAMS = params['score_params']
        self.HIT_PARAMS = params['hit_params']
        self.REQUIRED_USABLE_FRACTION = float(params['required_usable_fraction'])
        self.SEQUENCE_FILE_KEY = params['sequence_file_key']
        if isinstance(self.SEQUENCE_FILE_KEY, str):
            self.SEQUENCE_FILE_KEY = [self.SEQUENCE_FILE_KEY]
        self.LEVELS_TO_COPY = params['levels_to_copy']
        self.LEVEL_TO_USE_FOR_SCORES = params['level_to_use_for_scores']
        self.SCORE_KEY_FOR_TABLE = params['score_key_for_table']



class Metainfo:
        
    def __init__(self, meta_file):
        self.meta_file = meta_file
        if not os.path.exists(self.meta_file):
            with open(self.meta_file, 'w') as f:
                json.dump({}, f, indent=4)
        with open(self.meta_file) as f:
            self.META_PARAMS = json.load(f)

    def clear_meta_json(self):
        print(f'clearing meta json file {self.meta_file}')
        with open(self.meta_file, 'w') as f:
            json.dump({}, f, indent=4)

    def add_id_to_fail_dict(self, ref_id, reason: str):
        self.META_PARAMS['failed_id_dict'][str(ref_id)] = reason
        if str(ref_id) in self.META_PARAMS['ref_id_dict']:
            self.META_PARAMS['ref_id_dict'].pop(str(ref_id))

    def add_note_for_id(self, ref_id, note: str):
        id_entry = self.META_PARAMS['notes'].setdefault(str(ref_id), [])
        id_entry.append(note)

    def save_meta_dict(self):
        with open(self.meta_file, 'w') as f:
            json.dump(self.META_PARAMS, f, indent=4)

    def jsons_to_analyze(self):
        return list(self.META_PARAMS['ref_id_dict'].values())


