import json
import multiprocessing
from functools import partial
from pathlib import Path

import local_conservation_scores.PE_score_fl_alns as pe
import local_env_variables.env_variables as env

GAP_FRAC_CUTOFF = 0.2
OVERWRITE = True
LVL_DIR = env.ROOT / "data/example_orthogroup_database/human_odb_groups"
OUTPUT_DIR_NAME = "alignment_conservation_scores/property_entropy"
N_CORES = multiprocessing.cpu_count()


def run_pe(
    reference_id: str,
    alignment_file: str | Path,
    output_dir: str | Path,
    gap_frac_cutoff=GAP_FRAC_CUTOFF,
    overwrite=OVERWRITE,
):
    output_file = Path(output_dir) / Path(alignment_file).with_suffix(".json").name
    pe.main(
        alignment_file,
        output_file,
        gap_frac_cutoff=gap_frac_cutoff,
        overwrite=overwrite,
        reference_id=reference_id,
    )


def main(
    aln_filemap_file: str | Path,
    output_dir_name=OUTPUT_DIR_NAME,
    gap_frac_cutoff=GAP_FRAC_CUTOFF,
    overwrite=OVERWRITE,
    n_cores=N_CORES,
    multiprocess=True,
):
    aln_dir = Path(aln_filemap_file).parent
    with open(aln_filemap_file) as f:
        filemap = json.load(f)
    output_dir = Path(aln_dir) / output_dir_name
    output_dir.mkdir(parents=True, exist_ok=True)

    run_pe_w_params = partial(
        run_pe, gap_frac_cutoff=gap_frac_cutoff, overwrite=overwrite
    )
    if multiprocess:
        p = multiprocessing.Pool(n_cores)
        f_args = [(k, v, output_dir) for k, v in filemap.items()]
        p.starmap(run_pe_w_params, f_args)
        p.close()
        p.join()
    else:
        for k, v in filemap.items():
            run_pe_w_params(k, v, output_dir)


if __name__ == "__main__":
    # get a list of the directories in the LVL_DIR
    filemap_files = list(LVL_DIR.glob("**/filemap.json"))
    for filemap_file in filemap_files:
        print(filemap_file)
    for filemap_file in filemap_files:
        main(
            filemap_file,
            output_dir_name=OUTPUT_DIR_NAME,
            gap_frac_cutoff=GAP_FRAC_CUTOFF,
            overwrite=OVERWRITE,
            n_cores=N_CORES,
            multiprocess=True,
        )
