import json
import os
import re
import sys
from pathlib import Path

import logomaker as lm
import matplotlib.pyplot as plt

# plt.style.use("custom_standard")
# plt.style.use("custom_small")
import numpy as np
import pandas as pd
import seaborn as sns
from Bio import Align, AlignIO, Seq, SeqIO

import slim_conservation_scoring.seqtools.pssms as pssms


def build_mosaic_z_score_plot(figsize: tuple[int, int] = (15, 5)):
    mos_vector = []
    mos_vector.append(["background"] + ["scores"] * 2)
    mos_vector.append(["background"] + ["logo"] * 2)
    fig, axd = plt.subplot_mosaic(mos_vector, figsize=figsize, layout="constrained")
    return fig, axd


def plot_score_bar_plot(ax, score_list, query_seq):
    ax.bar(
        list(range(len(score_list))),
        score_list,
    )
    ax = _format_bar_plot(ax, query_seq)
    return ax


def _format_bar_plot(ax, xlabel_sequence: str):
    """format bar plot"""
    _ = ax.set_xticks(
        list(range(len(xlabel_sequence))),
        labels=list(xlabel_sequence),
    )
    ax.set_xlim(-0.5, len(xlabel_sequence) - 0.5)
    return ax


def plot_logo(ax, str_list, tick_label_str):
    counts = pssms.alignment_2_counts(str_list, show_plot=False, heatmap=False)
    lm.Logo(counts, color_scheme="chemistry", ax=ax)
    ax.set_ylim(0, len(str_list))
    _ = ax.set_xticks(
        list(range(len(str_list[0]))),
        labels=list(tick_label_str),
    )
    return ax


def format_logo_xticks_with_str(ax, tick_label_str):
    _ = ax.set_xticks(
        list(range(len(tick_label_str))),
        labels=list(tick_label_str),
    )
    return ax
