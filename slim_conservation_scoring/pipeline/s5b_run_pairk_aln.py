from pathlib import Path

import slim_conservation_scoring.pipeline.group_conservation_objects as group_tools
from slim_conservation_scoring.seqtools import general_utils as tools


from slim_conservation_scoring.conservation_scores import PairKmerAlnMethods

PAIRKMERALNFUNCS = PairKmerAlnMethods()


def flanked_pair_kmer_aln(
    og: group_tools.ConserGene,
    levels: list[str],
    score_key: str,
    function_name: str,
    function_params: dict,
    lflank: int,
    rflank: int,
    output_folder: str | Path | None = None,
):
    kmer_aln_function = PAIRKMERALNFUNCS.__getitem__(function_name)
    if output_folder is None:
        output_folder = Path(og.info_dict["analysis_folder"])
    output_folder = Path(output_folder)
    output_folder.mkdir(exist_ok=True, parents=True)
    # get idr sequence
    query_idr = og.query_sequence[og.idr_start : og.idr_end + 1]
    # get hit positions relative to the idr
    hit_st_idr = og.hit_start_position - og.idr_start
    hit_end_idr = og.hit_end_position - og.idr_start
    flanked_hit_st_idr, flanked_hit_end_idr, flanked_hit = tools.pad_hit(
        query_idr, hit_st_idr, hit_end_idr, lflank, rflank
    )
    orig_hit_st_in_flanked_hit = hit_st_idr - flanked_hit_st_idr
    orig_hit_end_in_flanked_hit = hit_end_idr - flanked_hit_st_idr
    k = len(flanked_hit)
    for level in levels:
        if level not in og.level_objects:
            continue
        lvlo = og.level_objects[level]
        aln_file = lvlo.alignment_file
        output_file = output_folder / f"{og.reference_index}-{level}-{score_key}.json"
        kmer_aln_function(
            input_alignment_file=aln_file,
            output_file=output_file,
            reference_id=og.query_gene_id,
            k=k,
            idr_aln_st=lvlo.idr_aln_start,
            idr_aln_end=lvlo.idr_aln_end,
            **function_params,
        )
        og.info_dict["orthogroups"][level]["conservation_scores"][f"{score_key}"] = {}
        score_dict = og.info_dict["orthogroups"][level]["conservation_scores"][
            f"{score_key}"
        ]
        score_dict["kmer_aln_file"] = str(output_file)
        score_dict["flanked_hit"] = flanked_hit
        score_dict["flanked_hit_start_position_in_idr"] = flanked_hit_st_idr
        score_dict["original_hit_st_in_flanked_hit"] = orig_hit_st_in_flanked_hit
        score_dict["original_hit_end_in_flanked_hit"] = orig_hit_end_in_flanked_hit
        score_dict["function_name"] = function_name
        score_dict["function_params"] = {
            k: v for k, v in function_params.items() if k != "mod"
        }
        score_dict["lflank"] = lflank
        score_dict["rflank"] = rflank
        og._overwrite_json()
    # gc.collect()


def pairwise_kmer_alignment_driver(
    json_file: str | Path,
    score_key: str,
    function_name: str,
    function_params: dict,
    levels: list[str] | None = None,
    lflank: int = 0,
    rflank: int = 0,
    score_output_folder: str | Path | None = None,
):
    og = group_tools.ConserGene(json_file)
    if hasattr(og, "critical_error"):
        return
    og.load_levels()
    if levels is None:
        levels = list(og.level_objects.keys())
    if hasattr(PAIRKMERALNFUNCS, function_name):
        flanked_pair_kmer_aln(
            og=og,
            levels=levels,
            score_key=score_key,
            function_name=function_name,
            function_params=function_params,
            lflank=lflank,
            rflank=rflank,
            output_folder=score_output_folder,
            # device=device,
            # threads=threads,
            # EsmMod=EsmMod,
        )
