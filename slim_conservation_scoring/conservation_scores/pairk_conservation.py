# %%
import concurrent.futures
import json
import multiprocessing
import time
from pathlib import Path
from typing import Callable

import slim_conservation_scoring.conservation_scores.tools.capra_singh_2007_scores as cs
from attrs import asdict, define, field, validators

# from alfpy import word_distance, word_pattern, word_vector
# from alfpy.utils import distmatrix
# from alfpy.utils import seqrecords as alf_seqrecords
import pairk


@define
class PairwiseScoreResults:
    hit_sequence: str
    hit_scores: list[float]
    hit_z_scores: list[float]
    background_scores: list[float]
    flank_hit_sequence: str | None = None
    flank_hit_scores: list[float] | None = None
    flank_hit_z_scores: list[float] | None = None


# def _kmer_distance_matrix(kmer_list, word_size=2):
#     alf_seq_records = alf_seqrecords.SeqRecords(id_list=kmer_list, seq_list=kmer_list)
#     p = word_pattern.create(alf_seq_records.seq_list, word_size=word_size)
#     # counts = word_vector.Counts(alf_seq_records.length_list, p)
#     freqs = word_vector.Freqs(alf_seq_records.length_list, p)
#     dist = word_distance.Distance(freqs, "google")
#     matrix = distmatrix.create(alf_seq_records.id_list, dist)
#     return matrix


# def _alfpy_query_matrix(refid: str, matrix) -> tuple[list[str], list[float]]:
#     """
#     get the row from the alfpy matrix that corresponds to the query gene
#     """
#     query_row = [c for c, i in enumerate(matrix.id_list) if refid in i][0]
#     query_row_distance = matrix.data[query_row]
#     query_row_similarity = 1 - query_row_distance
#     return matrix.id_list, query_row_similarity


# def _filter_similar_kmers(
#     kmer_list: list[str], hit_seq: str, similarity_threshold: float = 0.5
# ):
#     distmat = _kmer_distance_matrix(kmer_list, word_size=2)
#     id_list, sim = _alfpy_query_matrix(hit_seq, distmat)
#     passing_kmers = [i for i, s in zip(id_list, sim) if s <= similarity_threshold]
#     if hit_seq not in passing_kmers:
#         passing_kmers.append(hit_seq)
#     return passing_kmers


def pairk_conservation_from_json(
    kmer_aln_json: str | Path,
    hit_position: int,
    columnwise_score_func: Callable = cs.property_entropy,
    bg_cutoff: int = 50,
    bg_kmer_cutoff: int = 10,
) -> PairwiseScoreResults:
    """Calculate conservation scores from the pairk alignment results using the `pairk.calculate_conservation` function. Only the scores from the hit k-mer are returned.

    Parameters
    ----------
    kmer_aln_json : str | Path
        The json file storing the pairk alignment results
    hit_position : int
        The position of the hit k-mer in the query IDR sequence
    columnwise_score_func : Callable, optional
        A function to calculate conservation scores in a columnwise manner, can
        be any function that takes a string representing a column of an MSA.
        by default it is the property_entropy function from Capra and Singh 2007,
        DOI: 10.1093/bioinformatics/btm270.
    bg_cutoff : int, optional
        the minimum number of background scores required to calculate the
        z-scores, by default 50
    bg_kmer_cutoff : int, optional
        the minimum number of background kmers required to calculate the
        z-scores, by default 10

    Returns
    -------
    PairwiseScoreResults
        Object containing pairk conservation results for the hit k-mer.
        Includes the hit k-mer sequence, the conservation scores, the conservation
        z-scores, and the background scores.

    Raises
    ------
    ValueError
        If there are fewer than `bg_kmer_cutoff` kmers to use for background scores
    ValueError
        If there are fewer than `bg_cutoff` background scores
    """
    kmer_aln_json = Path(kmer_aln_json)
    pkaln = pairk.PairkAln.from_file(kmer_aln_json)
    if len(pkaln.orthokmer_matrix) < bg_kmer_cutoff:
        raise ValueError(
            f"not enough kmers in the matrices: {len(pkaln.orthokmer_matrix)}. Need at least {bg_kmer_cutoff}"
        )
    pkcons = pairk.calculate_conservation(pkaln, score_func=columnwise_score_func)
    if pkcons.n_bg_scores < bg_cutoff:
        raise ValueError(f"not enough background scores: {pkcons.n_bg_scores}")
    # pkcons.orthokmer_arr[hit_position, 0]
    assert (
        pkaln.orthokmer_matrix.loc[hit_position, "query_kmer"]
        == pkcons.orthokmer_arr[hit_position, 0]
    ), "hit kmer sequence mismatch between pairkaln and pairkcons"
    scores = PairwiseScoreResults(
        hit_sequence=pkcons.orthokmer_arr[hit_position, 0],
        hit_scores=list(pkcons.score_arr[hit_position, :]),
        hit_z_scores=list(pkcons.z_score_arr[hit_position, :]),
        background_scores=list(pkcons.bg_scores),
    )
    return scores
