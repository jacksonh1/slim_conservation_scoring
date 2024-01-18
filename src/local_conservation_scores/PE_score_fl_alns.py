#!/usr/bin/env python

# %%
import argparse
import copy
import json
import multiprocessing
from argparse import RawTextHelpFormatter
from pathlib import Path

import pandas as pd
from Bio import AlignIO, Seq, SeqIO
from pyprojroot import here

from local_conservation_scores import \
    conservation_scores_2007_capra_singh as cs
from local_conservation_scores import \
    conservation_scoring_tools as cons_tools


def score_alignment(alignment: AlignIO.MultipleSeqAlignment):
    aln = copy.deepcopy(alignment)
    scores = []
    for i in range(aln.get_alignment_length()):
        col = aln[:, i]
        scores.append(cs.property_entropy(col))
    return scores


def mask_alignment(alignment: AlignIO.MultipleSeqAlignment, reference_seq_str: str, gap_frac_cutoff: float = 0.2):
    gap_mask, score_mask = cons_tools.make_score_mask(
        alignment,
        reference_seq_str,
        gap_frac_cutoff=gap_frac_cutoff,
    )
    gap_mask = [bool(i) for i in gap_mask]
    score_mask = [bool(i) for i in score_mask]
    return gap_mask, score_mask


def main(
    input_alignment_file,
    output_file,
    gap_frac_cutoff,
    overwrite=False,
    reference_id=None,
):
    with open(input_alignment_file, 'r') as f:
        alignment = AlignIO.read(f, 'fasta')

    # check if output file exists
    if Path(output_file).exists() and not overwrite:
        print(f'{output_file} exists and overwrite is False')
        print('exiting...')
        return

    alignment_seqrecord_dict = {seqrecord.id: seqrecord for seqrecord in alignment}

    if reference_id is None:
        reference_id = alignment[0].id
    reference_seqrecord = alignment_seqrecord_dict[reference_id]
    reference_seq_str = str(reference_seqrecord.seq)
    score_dict = {}
    score_dict["gap_mask"], score_dict["score_mask"] = mask_alignment(
        alignment,
        reference_seq_str,
        gap_frac_cutoff=gap_frac_cutoff,
    )
    score_dict["scores"] = score_alignment(
        alignment,
    )

    with open(output_file, 'w') as f:
        json.dump(score_dict, f, indent=4)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='calculate asymmetric sum of pairs conservation scores (valdar) for an input alignment file',
        formatter_class=RawTextHelpFormatter
    )
    parser.add_argument(
        '-i',
        '--input',
        type=str,
        metavar='<file>',
        required=True,
        help='input alignment file (fasta format)'
    )
    parser.add_argument(
        '-o',
        '--output_file',
        type=str,
        metavar='<file>',
        required=True,
        help='output file (json)'
    )
    parser.add_argument(
        '--overwrite',
        action='store_true',
        help='overwrite existing results'
    )
    parser.add_argument(
        '-g',
        '--gap_frac_cutoff',
        type=float,
        metavar='<float>',
        default=0.2,
        help='''fraction of gaps allowed in a column. If column has >gap_frac_cutoff, it is masked in the gap mask, default=0.2'''
    )
    parser.add_argument(
        '-r',
        '--reference_id',
        type=str,
        metavar='<str>',
        default=None,
        help='''reference id to calculate score mask with (gaps in this sequence will be masked). If not provided, the first sequence in the alignment will be used'''
    )
    args = parser.parse_args()
    main(
        args.input,
        args.output_file,
        args.gap_frac_cutoff,
        overwrite=args.overwrite,
        reference_id=args.reference_id,
    )
