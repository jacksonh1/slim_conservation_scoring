
from pathlib import Path

import local_conservation_analysis_pipeline.group_conservation_objects as group_tools
from local_conservation_scores import conservation_score_methods

SCORES = conservation_score_methods()

def compute_conservation_scores(
    json_file: str|Path,
    score_key: str,
    score_params: dict,
):
    og = group_tools.ConserGene(json_file)
    og.load_levels()
    score_function = SCORES.__getitem__(score_key)

    output_folder = Path(og.info_dict["analysis_folder"])
    for level in og.level_objects:
        lvlo = og.level_objects[level]
        aln_file = lvlo.alignment_file
        output_file = output_folder / f"{aln_file.stem}-{score_key}.json"
        score_function(
            input_alignment_file=aln_file,
            output_file=output_file,
            reference_id=og.query_gene_id,
            **score_params,
        )
        
        og.info_dict["orthogroups"][level]['conservation_scores'][f"{score_key}"] = str(output_file)
        og._overwrite_json()
    
