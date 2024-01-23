
from pathlib import Path

import local_conservation_analysis_pipeline.group_conservation_objects as group_tools
import local_env_variables.env_variables as env


def compute_conservation_scores(
    json_file: str|Path,
    score_key: str,
    score_params: dict,
):
    og = group_tools.ConserGene(json_file)
    og.load_levels()

    output_folder = Path(og.info_dict["analysis_folder"])
    for level in og.level_objects:
        lvlo = og.level_objects[level]
        aln_file = lvlo.alignment_file
        output_file = output_folder / f"{aln_file.stem}-{score_key}.json"
        env.CONSERVATION_SCORE_METHODS[score_key](
            input_alignment_file=aln_file,
            output_file=output_file,
            reference_id=og.query_gene_id,
            **score_params,
        )
        
        og.info_dict["orthogroups"][level]['conservation_scores'][f"{score_key}_file"] = str(output_file)
        og._overwrite_json()
    




