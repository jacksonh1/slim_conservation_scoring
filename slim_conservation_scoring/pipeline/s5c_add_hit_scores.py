from pathlib import Path
import numpy as np

import slim_conservation_scoring.pipeline.group_conservation_objects as group_tools
from slim_conservation_scoring.conservation_scores import PairKmerConservationMethods
from slim_conservation_scoring.conservation_scores import ColumnwiseScoreMethods
from slim_conservation_scoring.conservation_scores import PairKmerAlnMethods

# import gc
from typing import Callable
from slim_conservation_scoring.conservation_scores.tools import (
    capra_singh_2007_scores as cs,
)
from attrs import asdict, define, field, validators
from slim_conservation_scoring.config import conservation_pipeline_parameters as conf


PAIR_KMER_CONS_FUNCS = PairKmerConservationMethods()
PAIR_KMER_ALN_FUNCS = PairKmerAlnMethods()
COLSCOREMETHODS = ColumnwiseScoreMethods()


@define
class PairwiseScoreResults:
    hit_sequence: str
    hit_scores: list[float]
    hit_z_scores: list[float]
    background_scores: list[float]
    flank_hit_sequence: str | None = None
    flank_hit_scores: list[float] | None = None
    flank_hit_z_scores: list[float] | None = None


def lvlo_2_pairwise_scores(
    lvlo: group_tools.ConserLevel,
    score_key: str,
    params: conf.PairKmerConservationParams,
):
    kmer_conservation_function = PAIR_KMER_CONS_FUNCS.__getitem__(
        params.kmer_conservation_function_name
    )
    col_function = COLSCOREMETHODS.__getitem__(params.columnwise_score_function_name)

    pairdict = lvlo.conservation_scores[score_key]
    flanked_hit_scores = kmer_conservation_function(
        pairdict["kmer_aln_file"],
        pairdict["flanked_hit_start_position_in_idr"],
        columnwise_score_func=col_function,
        bg_cutoff=params.bg_cutoff,
        bg_kmer_cutoff=params.bg_kmer_cutoff,
    )
    hit_slice = slice(
        pairdict["original_hit_st_in_flanked_hit"],
        pairdict["original_hit_end_in_flanked_hit"] + 1,
    )
    scores = PairwiseScoreResults(
        hit_sequence=flanked_hit_scores.hit_sequence[hit_slice],
        hit_scores=flanked_hit_scores.hit_scores[hit_slice],
        hit_z_scores=flanked_hit_scores.hit_z_scores[hit_slice],
        flank_hit_sequence=flanked_hit_scores.hit_sequence,
        flank_hit_scores=flanked_hit_scores.hit_scores,
        flank_hit_z_scores=flanked_hit_scores.hit_z_scores,
        background_scores=flanked_hit_scores.background_scores,
    )
    return scores


def pairwise_scores(
    og: group_tools.ConserGene,
    levels: list[str],
    score_key: str,
    params: conf.PairKmerConservationParams,
):
    for level in levels:
        if level not in og.level_objects:
            continue
        lvlo = og.level_objects[level]
        if score_key not in lvlo.conservation_scores:
            raise ValueError(f"score_key {score_key} not in lvlo.conservation_scores")
        score_dict = og.info_dict["orthogroups"][level]["conservation_scores"][
            f"{score_key}"
        ]
        try:
            scores = lvlo_2_pairwise_scores(
                lvlo=lvlo,
                score_key=score_key,
                params=params,
            )
        except ValueError as e:
            score_dict["score_error"] = str(e)
            continue
        score_dict["flanked_hit_sequence"] = scores.flank_hit_sequence
        score_dict["flanked_hit_scores"] = scores.flank_hit_scores
        score_dict["flanked_hit_z_scores"] = scores.flank_hit_z_scores
        score_dict["hit_sequence"] = scores.hit_sequence
        score_dict["hit_scores"] = scores.hit_scores
        score_dict["hit_z_scores"] = scores.hit_z_scores
        score_dict["bg_std"] = np.std(scores.background_scores)
        score_dict["pairk_conservation_params"] = asdict(params)
        og._overwrite_json()
    # gc.collect()


def compute_hit_conservation_scores(
    json_file: str | Path,
    pairk_method: conf.PairKAlnMethod | conf.PairKEmbeddingAlnMethod,
    params: conf.PairKmerConservationParams,
):
    og = group_tools.ConserGene(json_file)
    if hasattr(og, "critical_error"):
        return
    og.load_levels()
    if pairk_method.level is not None:
        levels = [pairk_method.level]
    else:
        levels = list(og.level_objects.keys())
    if hasattr(PAIR_KMER_ALN_FUNCS, pairk_method.function_name):
        pairwise_scores(
            og=og,
            levels=levels,
            score_key=pairk_method.score_key,
            params=params,
        )
