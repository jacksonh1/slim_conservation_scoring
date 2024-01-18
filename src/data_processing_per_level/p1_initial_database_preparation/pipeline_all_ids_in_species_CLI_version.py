import argparse
import multiprocessing
import shutil
import traceback
from pathlib import Path

import local_config.conf as conf
import local_orthoDB_group_tools.sql_queries as sql_queries

import local_scripts.create_alnmap_no_levels as create_filemap
import local_scripts.odb_group_pipeline as pipeline

SPECIES_ID = "9606_0"
N_CORES = multiprocessing.cpu_count() - 2
OG_LEVELS = ["Eukaryota", "Mammalia", "Metazoa", "Tetrapoda", "Vertebrata"]


def run_pipeline(
    config: conf.PipelineParams, query_odb_gene_id: str
):
    '''
    run the pipeline for a single odb_gene_id for multiple og_levels
    '''
    try:
        pipeline.main_pipeline(config, odb_gene_id=query_odb_gene_id)
    except ValueError as err:
        traceback.print_exc()
        # logger.error(f"{query_geneid} - {og_level} - {err}")
        print(f"{query_odb_gene_id} - {og_level} - {err}")


def main(
    config: conf.PipelineParams,
    multiprocess=True,
    species_id=SPECIES_ID,
    n_cores=N_CORES,
    overwrite=False,
):
    odbgeneid_list = sql_queries.get_all_odb_gene_ids_from_species_id(species_id)
    if Path(config.main_output_folder).exists():
        if overwrite:
            shutil.rmtree(config.main_output_folder)
        else:
            raise FileExistsError(
                f"main_output_folder already exists: {config.main_output_folder}. Use -o flag to overwrite"
            )
    if multiprocess:
        p = multiprocessing.Pool(n_cores)
        f_args = [(config, i) for i in odbgeneid_list]
        p.starmap(run_pipeline, f_args)
        p.close()
        p.join()
    else:
        for i in odbgeneid_list:
            run_pipeline(config, i)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="run the pipeline for all genes in an organism. Create a new folder for each og_level in the main_output_folder",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-c",
        "--config",
        type=str,
        metavar="<file>",
        default=None,
        help="""path to config file""",
    )
    parser.add_argument(
        "-n",
        "--n_cores",
        type=int,
        metavar="<int>",
        default=N_CORES,
        help=f"""number of cores to use""",
    )
    parser.add_argument(
        "-s",
        "--species_id",
        type=str,
        metavar="<str>",
        default=SPECIES_ID,
        help=f"""species id to use""",
    )
    parser.add_argument(
        "-o",
        "--overwrite",
        action="store_true",
        help="""if flag is provided and the main_output_folder exists, it will be removed and overwritten by the new files. Otherwise, an error will be raised if the folder exists""",
    )
    args = parser.parse_args()
    config = pipeline.load_config(args.config)
    parent_folder = Path(config.main_output_folder)
    for og_level in OG_LEVELS:
        config.main_output_folder = str(parent_folder / og_level)
        config.og_select_params.OG_level_name = og_level
        main(
            config,
            multiprocess=True,
            species_id=args.species_id,
            n_cores=args.n_cores,
            overwrite=args.overwrite,
        )
        create_filemap.create_filemap(
            main_output_folder=config.main_output_folder,
            output_file=Path(config.main_output_folder) / "filemap.json",
        )
