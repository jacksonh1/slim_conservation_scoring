# %%
import json
import multiprocessing
from collections import defaultdict
from pathlib import Path

# from loguru import logger
from pyprojroot import here

import local_conservation_analysis_pipeline.searching_for_hit_sequence.param_and_log_tools as param_and_log_tools
from local_orthoDB_analysis_tools_v2 import hit_sequence_tools as hit_tools
from local_seqtools import general_utils as tools


def prep_og_hit_class_at_level(
        og_json_file,
        hit_sequence_in,
        longest_common_substring = False,
        LCS_min_length=10,
        target_hit_length=36
    ):
    ogseqs_o = hit_tools.hit_sequence_class(
        og_json_file
    )
    if longest_common_substring:
        hit_sequence = tools.longest_common_substring(
            hit_sequence_in,
            ogseqs_o.query_sequence_str,
        )
        if len(hit_sequence) < LCS_min_length:
            ogseqs_o.found_hit = False
            return ogseqs_o
    else:
        hit_sequence = hit_sequence_in
    ogseqs_o.find_query_hit_sequence_in_alignment(
        hit_sequence,
        padding_target_length=target_hit_length
    )
    if not ogseqs_o.found_hit:
        return ogseqs_o
    ogseqs_o.hit_idr_characterization(cider=False)
    ogseqs_o._slice_alignment()
    return ogseqs_o


def add_info_2_analysis_dict(og_info, ogseqs_o):
    og_info['found hit'] = ogseqs_o.found_hit
    og_info['query_sequence_id_str'] = ogseqs_o.query_sequence_id_str
    og_info['query_sequence_str'] = ogseqs_o.query_sequence_str
    if not ogseqs_o.found_hit:
        return og_info
    og_info['hit_sequence'] = ogseqs_o.hit_sequence
    og_info['hit_start_position'] = ogseqs_o.hit_start_position
    og_info['hit_end_position'] = ogseqs_o.hit_end_position
    og_info['hit_in_idr'] = ogseqs_o.hit_in_idr
    if ogseqs_o.hit_in_idr:
        og_info['idr_start'] = ogseqs_o.idr_start
        og_info['idr_end'] = ogseqs_o.idr_end
    return og_info


def driver(og_info_json, params_in, logger, root=here()):
    params = params_in.copy()
    og_info_json = og_info_json
    with open(og_info_json, 'r') as f:
        og_info = json.load(f)
    logger.debug(f"generating json files with hit information for {og_info['query_gene_id']}-{og_info['UniprotID']}-{og_info['orig_hit_sequence']}")
    analysis_folder = str(og_info_json.resolve().relative_to(root).parent)
    og_info['analysis_folder'] = analysis_folder
    new_levels_list = []
    og_info['level_info'] = defaultdict(dict)
    for level in og_info['levels']:
        level_info_dict = {}
        ogseqs_o = prep_og_hit_class_at_level(
            og_info['json_files_per_level'][level],
            og_info['orig_hit_sequence'],
            **params
        )
        og_info = add_info_2_analysis_dict(og_info, ogseqs_o)
        if len(ogseqs_o.alignment_orthologs_ldo_clustered) < 2:
            logger.warning(f"insufficient orthologs for {og_info['reference_index']}-{og_info['query_gene_id']}-{og_info['UniprotID']} at {level} level")
            continue
        if not ogseqs_o.found_hit:
            logger.warning(f"hit sequence {og_info['orig_hit_sequence']} not found for {og_info['reference_index']}-{og_info['query_gene_id']}-{og_info['UniprotID']}. Or LCS too short")
            continue
        new_levels_list.append(level)
        level_info_dict=ogseqs_o.add_attributes2dict(level_info_dict)
        level_info_dict['UniprotID'] = og_info['UniprotID']
        level_info_dict['name'] = og_info['name']
        level_info_dict['reference_index'] = og_info['reference_index']
        level_info_dict['analysis_folder'] = analysis_folder
        og_info['level_info'][level] = level_info_dict
        # alignment_info_json_name = og_info_json.parent / f"{level}_aln_hit_info_.json"
        # with open(alignment_info_json_name, 'w') as f:
            # json.dump(level_info_dict, f, indent=4)
    og_info['levels'] = new_levels_list
    with open(og_info_json, 'w') as f:
        json.dump(og_info, f, indent=4)

def main(
    score_params: dict,
    logger,
    meta_info: param_and_log_tools.Metainfo,
    root: Path,
    multiprocess=True
):
    og_info_jsons = meta_info.jsons_to_analyze()
    og_info_jsons = [root / i for i in og_info_jsons]
    if multiprocess:
        p = multiprocessing.Pool(
            multiprocessing.cpu_count() - 1
        )
        f_args = [(i, score_params, logger, root) for i in og_info_jsons]
        p.starmap(driver, f_args)
        # might be better to use imap_unordered or something that will close each process as it finishes
        # p.map(driver, og_info_jsons)
        p.close()
        p.join()
    else:
        # driver([i for i in og_info_jsons if '9606_0:0014fb' in i][0], PARAMS)
        for i in og_info_jsons:
            print(i)
            driver(i, score_params, logger, root)



# %%
