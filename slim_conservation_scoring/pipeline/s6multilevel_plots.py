# %%
import json
from pathlib import Path

import matplotlib.pyplot as plt
from Bio import AlignIO, SeqIO

import slim_conservation_scoring.pipeline.group_conservation_objects as group_tools
import slim_conservation_scoring.conservation_scores.tools.score_plots as score_plots
from slim_conservation_scoring.config import conservation_pipeline_parameters as conf
import pairk
import slim_conservation_scoring.env_variables.env_variables as env


from slim_conservation_scoring.conservation_scores import MSAScoreMethods
from slim_conservation_scoring.conservation_scores import PairKmerAlnMethods
from slim_conservation_scoring.conservation_scores import ColumnwiseScoreMethods

MSASCORES = MSAScoreMethods()
PAIRKALNFUNCS = PairKmerAlnMethods()
COLSCOREMETHODS = ColumnwiseScoreMethods()


# %%


def make_msa_multi_level_plot(json_file, config: conf.PipelineParameters):
    og = group_tools.ConserGene(json_file)
    og.load_aln_scores(
        config.multilevel_plot_params.score_key,
        num_bg_scores_cutoff=config.multilevel_plot_params.num_bg_scores_cutoff,
    )
    output_folder = Path(og.info_dict["analysis_folder"])
    multi_plot_filename = score_plots.consergene_2_multilevel_plot(
        og,
        big_plot_output_folder=output_folder,
        score_name=f"{config.multilevel_plot_params.score_key}",
        score_type=config.multilevel_plot_params.score_type,
        bar_ylim=[-0.5, 2.5],
        strip_gaps=config.multilevel_plot_params.strip_gaps,
    )
    og.add_item_to_json(
        f"multilevel_plot_file-{config.multilevel_plot_params.score_key}",
        str(multi_plot_filename),
        save_json=True,
    )
    plt.close("all")


def make_pairk_multilevel_plot(json_file, config: conf.PipelineParameters):
    col_function = COLSCOREMETHODS.__getitem__(
        config.pairk_conservation_params.columnwise_score_function_name
    )
    og = group_tools.ConserGene(json_file)
    og.load_levels()
    levels = list(og.info_dict["orthogroups"].keys())
    if all([x in env.PHYLOGENY_LVL_ORDERING for x in levels]):
        levels = sorted(
            levels,
            key=lambda x: env.PHYLOGENY_LVL_ORDERING.index(x),
        )
    plt.rcParams["font.size"] = 10
    fig, axd = score_plots.build_og_level_screen_mosaic_z_score(levels)
    output_folder = Path(og.info_dict["analysis_folder"])
    for level in levels:
        if level not in og.levels_passing_filters:
            message = f"{og.query_gene_id}: not enough sequences for {level}"
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
        lvlo = og.level_objects[level]
        if config.multilevel_plot_params.score_key not in lvlo.conservation_scores:
            message = f"{og.query_gene_id}: no scores were calculated at the {level} level for {config.multilevel_plot_params.score_key}"
            axd[f"scores-{level}"].text(
                0.5,
                0.5,
                message,
                horizontalalignment="center",
                verticalalignment="center",
                transform=axd[f"scores-{level}"].transAxes,
                fontsize=11,
            )
            continue
        kmer_aln_json = lvlo.conservation_scores[
            config.multilevel_plot_params.score_key
        ]["kmer_aln_file"]
        fl_hit_start = lvlo.conservation_scores[
            config.multilevel_plot_params.score_key
        ]["flanked_hit_start_position_in_idr"]
        pkaln = pairk.PairkAln.from_file(kmer_aln_json)
        pkcons = pairk.calculate_conservation(pkaln, score_func=col_function)
        if level not in og.levels_passing_filters:
            message = f"{og.query_gene_id}: not enough sequences"
            axd[f"scores-{level}"].text(
                0.5,
                0.5,
                message,
                horizontalalignment="center",
                verticalalignment="center",
                transform=axd[f"scores-{level}"].transAxes,
                fontsize=11,
            )
            continue
        if pkcons.n_bg_scores < config.multilevel_plot_params.num_bg_scores_cutoff:
            message = f"{og.query_gene_id}: not enough background scores for z-score ({pkcons.n_bg_scores} < cutoff ({config.multilevel_plot_params.num_bg_scores_cutoff}))"
            axd[f"scores-{level}"].text(
                0.5,
                0.5,
                message,
                horizontalalignment="center",
                verticalalignment="center",
                transform=axd[f"scores-{level}"].transAxes,
                fontsize=11,
            )
            pkcons.plot_background_distribution(axd[f"bg_dist-{level}"])
            continue
        pkcons.plot_background_distribution(axd[f"bg_dist-{level}"])
        pkcons.plot_score_barplot(
            position=fl_hit_start,
            score_type=config.multilevel_plot_params.score_type,
            ax=axd[f"scores-{level}"],
        )
        pkcons.plot_sequence_logo(fl_hit_start, axd[f"logo-{level}"])
        plot_title = f"{og.query_gene_id} - {level} - {pkcons.orthokmer_arr.shape[1]} sequences - z-score (bg of {pkcons.n_bg_scores} scores from {pkcons.n_bg_kmers} k-mers)"  # type: ignore
        axd[f"scores-{level}"].set_title(plot_title, fontsize=11)
        axd[f"bg_dist-{level}"].set_title(f"{level} - {og.query_gene_id}", fontsize=11)
        axd[f"bg_dist-{level}"].set_xlabel("")
    output_file = (
        output_folder
        / f"{og.reference_index}-{og.query_gene_id.replace(':','')}-{config.multilevel_plot_params.score_key}_multilevel_plot.png"
    )
    # plt.tight_layout(h_pad=0.5, w_pad=0.5)
    # plt.subplots_adjust(wspace=0.4, hspace=0.4)
    # fig.savefig(output_file, dpi=300)
    fig.savefig(output_file, bbox_inches="tight", dpi=300)
    og.add_item_to_json(
        f"multilevel_plot_file-{config.multilevel_plot_params.score_key}",
        str(output_file),
        save_json=True,
    )
    plt.close("all")


def multi_level_plot_driver(json_file, config: conf.PipelineParameters):
    scorekey_dict = config.get_score_key_dict()
    if config.multilevel_plot_params.score_key is None:
        return
    if config.multilevel_plot_params.score_key not in scorekey_dict:
        print(
            f"multilevel plot score_key {config.multilevel_plot_params.score_key} was not found in the config file"
        )
        print(f"Available score keys are: {scorekey_dict.keys()}")
        print("Skipping multilevel plot")
        return
    scoremethod = scorekey_dict[config.multilevel_plot_params.score_key]
    if hasattr(PAIRKALNFUNCS, scoremethod.function_name):
        make_pairk_multilevel_plot(json_file, config)
    elif hasattr(MSASCORES, scoremethod.function_name):
        make_msa_multi_level_plot(
            json_file,
            config,
        )


# lvlo.load_scores('property_entropy')

# dir(lvlo)

# group_tools.conservation_level_score(og.info_dict['orthogroups'][level])

# %%
