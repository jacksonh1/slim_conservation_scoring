import json
import os
import re
import sys
from pathlib import Path

import logomaker as lm
import matplotlib.pyplot as plt
import seaborn as sns

import local_seqtools.pssms as pssms
from local_orthoDB_analysis_tools_v2 import group_tools_v2 as group_tools

plt.style.use("custom_standard")
plt.style.use("custom_small")
import numpy as np
import pandas as pd
from Bio import Align, AlignIO, Seq, SeqIO


def format_bar_plot(ax, xlabel_sequence: str, bar_ylim=[0, 1], labelsize=16):
    """format bar plot"""
    _ = ax.set_xticks(
        list(range(len(xlabel_sequence))),
        labels=list(xlabel_sequence),
    )
    ax.set_xlim(-0.5, len(xlabel_sequence) - 0.5)
    ax.tick_params(axis="x", which="major", labelsize=labelsize)
    if bar_ylim is not None:
        ax.set_ylim(bar_ylim)
    return ax


def get_non_gap_indexes(hit_seq_in_alignment_str):
    """get list of nongap positions in `hit_seq_in_alignment_str`"""
    return [c for c, i in enumerate(hit_seq_in_alignment_str) if i != "-"]


def strip_gaps_from_slice(alignment_slice, query_seq_in_alignment_str):
    """strip gaps from alignment slice"""
    new_sequences = [""] * (len(alignment_slice))
    non_gap_indices = get_non_gap_indexes(query_seq_in_alignment_str)
    new_query_slice = "".join(
        [query_seq_in_alignment_str[i] for i in non_gap_indices]
    )
    for c in non_gap_indices:
        for j in range(len(new_sequences)):
            new_sequences[j] += alignment_slice[j, c]
    return new_sequences, new_query_slice, non_gap_indices


def plot_alignment_slice_conservation2query(
    query_seq_in_alignment_str: str,
    alignment: AlignIO.MultipleSeqAlignment,
    score_list: list,
    score_mask: list = None,
    slice_coordinates: list = None,
    strip_gaps: bool = True,
    bar_ylim=[0, 1],
    axes: list[plt.Axes] = None,
    figsize=(20, 5),
    **tick_params,
):
    """plot alignment slice conservation to query"""
    assert len(query_seq_in_alignment_str) == len(score_list)
    if axes is None:
        fig, (ax1, ax2) = plt.subplots(figsize=figsize, nrows=2)
    else:
        ax1, ax2 = axes
    if score_mask is not None:
        assert len(score_list) == len(score_mask)
        score_list = np.array(score_list)
        score_list[[not i for i in score_mask]] = 0
        score_list = list(score_list)
    if slice_coordinates is not None:
        alignment = alignment[:, slice_coordinates[0] : slice_coordinates[1]+1]
        query_seq_in_alignment_str = query_seq_in_alignment_str[
            slice_coordinates[0] : slice_coordinates[1]+1
        ]
        score_list = score_list[slice_coordinates[0] : slice_coordinates[1]+1]
        score_mask = score_mask[slice_coordinates[0] : slice_coordinates[1]+1]
    if strip_gaps:
        (
            alignment_str_list,
            query_seq_in_alignment_str,
            non_gap_indices,
        ) = strip_gaps_from_slice(
            alignment, query_seq_in_alignment_str
        )
        score_list = [score_list[i] for i in non_gap_indices]
    else:
        alignment_str_list = [str(i.seq) for i in alignment]
    if bar_ylim == "auto":
        limit = max(abs(min(score_list)), abs(max(score_list)))
        bar_ylim = [-limit, limit]

    ax1.bar(
        list(range(len(query_seq_in_alignment_str))),
        score_list,
    )
    ax1 = format_bar_plot(
        ax1, query_seq_in_alignment_str, bar_ylim=bar_ylim, **tick_params
    )

    counts = pssms.alignment_2_counts(
        alignment_str_list, show_plot=False, heatmap=False
    )
    lm.Logo(counts, color_scheme="chemistry", ax=ax2)
    ax2.set_ylim(0, len(alignment_str_list))
    return ax1, ax2, counts


def bg_score_distro_plot(score_list: list, ax):
    sns.histplot(score_list, bins=20, ax=ax)
    ax.set_xlim([0, 1])
    # add a vertical line at the mean
    ax.axvline(np.mean(score_list), color="k", linewidth=1)
    # add a vertical line at one standard deviation above and below the mean
    ax.axvline(
        np.mean(score_list) + np.std(score_list),
        color="k",
        linestyle="dashed",
        linewidth=1,
    )
    ax.axvline(
        np.mean(score_list) - np.std(score_list),
        color="k",
        linestyle="dashed",
        linewidth=1,
    )

# ==============================================================================
# // z-score plots
# ==============================================================================
def build_og_level_screen_mosaic_z_score(levels):
    mos_vector = []
    for level in levels:
        mos_vector.append([f"bg_dist-{level}"] + [f"scores-{level}"] * 3)
        mos_vector.append([f"bg_dist-{level}"] + [f"logo-{level}"] * 3)
    fig, axd = plt.subplot_mosaic(mos_vector, figsize=(10, 10), layout="constrained")
    plt.tight_layout()
    return fig, axd


def add_z_score_plots_to_mosaic(
    axd, level, og_scores_o: group_tools.level_score_class, bar_ylim=[-2.5, 2.5], highlight_positions=None, **plot_kwargs
):
    ax1, ax2, counts = plot_alignment_slice_conservation2query(
        query_seq_in_alignment_str=og_scores_o.query_aligned_sequence_str,
        alignment=og_scores_o.alignment_orthologs_ldo_clustered,
        score_list=og_scores_o.z_score_dict["z_scores"],
        score_mask=og_scores_o.score_mask,
        slice_coordinates=[
            og_scores_o.hit_alignment_start_position,
            og_scores_o.hit_alignment_end_position,
        ],
        bar_ylim=bar_ylim,
        axes=[axd[f"scores-{level}"], axd[f"logo-{level}"]],
        **plot_kwargs,
    )
    if highlight_positions is not None:
        for position in highlight_positions:
            for ax in [ax1, ax2]:
                ax.axvspan(
                    position+0.5,
                    position-0.5,
                    color="red",
                    alpha=0.3,
                )


def big_z_score_plot(og_folder: group_tools.ortholog_group_info, **kwargs):
    fig, axd = build_og_level_screen_mosaic_z_score(og_folder.levels)
    for level in og_folder.levels:
        og_scores_o = og_folder.level_score_objects[level]
        if og_scores_o.z_score_dict is None:
            print(f"no bg scores for {og_folder.name}, {level}, skipping")
            message = f"{og_folder.name} ({og_folder.UniprotID}) - {level} - {len(og_scores_o.alignment_orthologs_ldo_clustered)} sequences - no bg scores"
            axd[f"scores-{level}"].text(
                0.5,
                0.5,
                message,
                horizontalalignment="center",
                verticalalignment="center",
                transform=axd[f"scores-{level}"].transAxes,
                fontsize=11,
            )
            # axd[f'scores-{level}'].set_title(, fontsize=11)
            continue
        add_z_score_plots_to_mosaic(axd, level, og_scores_o, **kwargs)
        bg_score_distro_plot(
            og_scores_o.z_score_dict["bg_scores"], axd[f"bg_dist-{level}"]
        )
        plot_title = f"{og_folder.name} ({og_folder.UniprotID}) - {level} - {len(og_scores_o.alignment_orthologs_ldo_clustered)} sequences - z-score (IDR bg using {len(og_scores_o.z_score_dict['bg_scores'])} residues)"
        axd[f"scores-{level}"].set_title(plot_title, fontsize=11)
        axd[f"bg_dist-{level}"].set_title(f"{level} - {og_folder.UniprotID}")
    return fig, axd


# ==============================================================================
# // just score plots (not z-scores)
# ==============================================================================
def build_og_level_screen_mosaic(levels):
    mos_vector = []
    for level in levels:
        mos_vector.append([f"scores-{level}"] * 3)
        mos_vector.append([f"logo-{level}"] * 3)
    fig, axd = plt.subplot_mosaic(mos_vector, figsize=(10, 10), layout="constrained")
    plt.tight_layout()
    return fig, axd


def add_score_plots_to_mosaic(
    axd, level, og_scores_o, bar_ylim=[0, 1], highlight_positions=None, **plot_kwargs
):
    ax1, ax2, counts = plot_alignment_slice_conservation2query(
        query_seq_in_alignment_str=og_scores_o.query_aligned_sequence_str,
        alignment=og_scores_o.alignment_orthologs_ldo_clustered,
        score_list=og_scores_o.scores,
        score_mask=og_scores_o.score_mask,
        slice_coordinates=[
            og_scores_o.hit_alignment_start_position,
            og_scores_o.hit_alignment_end_position,
        ],
        bar_ylim=bar_ylim,
        axes=[axd[f"scores-{level}"], axd[f"logo-{level}"]],
        **plot_kwargs,
    )
    if highlight_positions is not None:
        for position in highlight_positions:
            for ax in [ax1, ax2]:
                ax.axvspan(
                    position+0.5,
                    position-0.5,
                    color="red",
                    alpha=0.3,
                )


def big_score_plot(og_folder: group_tools.ortholog_group_info, **kwargs):
    fig, axd = build_og_level_screen_mosaic(og_folder.levels)
    for level in og_folder.levels:
        og_scores_o = og_folder.level_score_objects[level]
        add_score_plots_to_mosaic(axd, level, og_scores_o, **kwargs)
        plot_title = f"{og_folder.name} ({og_folder.UniprotID}) - {level} - {len(og_scores_o.alignment_orthologs_ldo_clustered)} sequences - NOT Z-SCORES. Hit not in IDR."
        axd[f"scores-{level}"].set_title(plot_title, fontsize=11)
    return fig, axd


# %%
# ==============================================================================
# // big plot driver
# ==============================================================================
def ogscreen_plot_from_og_folder_class(
    og_folder: group_tools.ortholog_group_info, save_big_plot=True, big_plot_output_folder=None, score_name=None, **kwargs
):
    if og_folder.hit_in_idr:
        fig, axd=big_z_score_plot(og_folder, **kwargs)
    else:
        fig, axd=big_score_plot(og_folder, **kwargs)
    if big_plot_output_folder is None:
        big_plot_output_folder = og_folder.analysis_folder
    big_plot_filename = (
        big_plot_output_folder
        / f"{og_folder.reference_index}-{og_folder.name}-{og_folder.UniprotID}-{og_folder.query_gene_id}_og_level_score_screen.png"
    )
    if save_big_plot:
        big_plot_output_folder = Path(big_plot_output_folder)
        big_plot_output_folder.mkdir(parents=True, exist_ok=True)
        if score_name is None:
            big_plot_filename = (
                big_plot_output_folder
                / f"{og_folder.reference_index}-{og_folder.name}-{og_folder.UniprotID}-{og_folder.query_gene_id}_og_level_score_screen.png"
            )
        else:
            big_plot_filename = (
                big_plot_output_folder
                / f"{og_folder.reference_index}-{og_folder.name}-{og_folder.UniprotID}-{og_folder.query_gene_id}-{score_name}_og_level_score_screen.png"
            )
        fig.savefig(big_plot_filename, bbox_inches="tight", dpi=300)
    return big_plot_filename
    # return fig, axd, big_plot_filename


# %%
# ==============================================================================
# // TITLE
# ==============================================================================
def build_mosaic_z_score_plot():
    mos_vector = []
    mos_vector.append(["bg_dist"] + ["scores"]*2)
    mos_vector.append(["bg_dist"] + ["logo"]*2)
    fig, axd = plt.subplot_mosaic(mos_vector, figsize=(15, 5), layout="constrained")
    plt.tight_layout()
    return fig, axd

def plot_score_bar_plot(ax, score_list, query_seq, mask=None):
    if mask is not None:
        score_list = np.array(score_list)
        score_list[[not i for i in mask]] = 0
        score_list = list(score_list)
    ax.bar(
        list(range(len(score_list))),
        score_list,
    )
    ax = _format_bar_plot(ax, query_seq, bar_ylim=False)
    return ax

def _format_bar_plot(ax, xlabel_sequence: str, bar_ylim=[0, 1], labelsize=16):
    """format bar plot"""
    _ = ax.set_xticks(
        list(range(len(xlabel_sequence))),
        labels=list(xlabel_sequence),
    )
    ax.set_xlim(-0.5, len(xlabel_sequence) - 0.5)
    ax.tick_params(axis="x", which="major", labelsize=labelsize)
    if bar_ylim:
        ax.set_ylim(bar_ylim)
    return ax

def plot_logo(ax, str_list, tick_label_str):
    counts = pssms.alignment_2_counts(
        str_list, show_plot=False, heatmap=False
    )
    lm.Logo(counts, color_scheme="chemistry", ax=ax)
    ax.set_ylim(0, len(str_list))
    _=ax.set_xticks(
        list(range(len(str_list[0]))),
        labels=list(tick_label_str),
    )
    return ax

def format_logo_xticks_with_str(ax, tick_label_str):
    _=ax.set_xticks(
        list(range(len(tick_label_str))),
        labels=list(tick_label_str),
    )
    return ax
