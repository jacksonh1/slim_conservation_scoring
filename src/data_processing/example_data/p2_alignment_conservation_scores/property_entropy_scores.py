import json
import multiprocessing
from functools import partial
from pathlib import Path

import local_conservation_scores.aln_property_entropy as pe

GAP_FRAC_CUTOFF = 0.2
OVERWRITE = True
database_dir = Path("../../../../data/example_orthogroup_database/human_odb_groups")
JSON_DIR = database_dir / "info_jsons"
OUTPUT_DIR = database_dir / "alignment_conservation_scores" / "aln_property_entropy"
N_CORES = multiprocessing.cpu_count()


def run_pe(data, output_dir, gap_frac_cutoff=GAP_FRAC_CUTOFF, overwrite=OVERWRITE):
    odb_gene_id = data["query_odb_gene_id"]
    alignment_file = Path(data["alignment_clustered_ldos_file"])
    output_file = output_dir / f"{alignment_file.stem}.json"
    if output_file.exists() and not overwrite:
        raise FileExistsError(f"{output_file} already exists and overwrite=False")
    if "conservation_scores" not in data:
        data["conservation_scores"] = {}
    try:
        pe.main(
            alignment_file,
            output_file,
            gap_frac_cutoff=gap_frac_cutoff,
            overwrite=overwrite,
            reference_id=odb_gene_id,
        )
    except ZeroDivisionError as err:
        data["conservation_scores"]["aln_property_entropy"] = f"error - {str(err)}"
        return data
    pe_dict = {}
    pe_dict["file"] = str(output_file.resolve())
    pe_dict["score_function_name"] = "aln_property_entropy"
    pe_dict["score_params"] = {
        "gap_frac_cutoff": gap_frac_cutoff,
        "overwrite": overwrite,
    }
    data["conservation_scores"]["aln_property_entropy"] = pe_dict
    return data


def run_pe_driver(
    json_file: str | Path,
    output_dir,
    gap_frac_cutoff=GAP_FRAC_CUTOFF,
    overwrite=OVERWRITE,
):
    with open(json_file) as f:
        data = json.load(f)
    try:
        data = run_pe(
            data, output_dir, gap_frac_cutoff=GAP_FRAC_CUTOFF, overwrite=OVERWRITE
        )
    except FileExistsError as ferr:
        print(str(ferr))
        return
    with open(json_file, "w") as f:
        json.dump(data, f, indent=4)


def main(
    json_dir=JSON_DIR,
    output_dir=OUTPUT_DIR,
    gap_frac_cutoff=GAP_FRAC_CUTOFF,
    overwrite=OVERWRITE,
    n_cores=N_CORES,
    multiprocess=True,
):
    output_dir.mkdir(parents=True, exist_ok=True)
    info_json_files = list(json_dir.glob("*.json"))

    run_pe_w_params = partial(
        run_pe_driver, gap_frac_cutoff=gap_frac_cutoff, overwrite=overwrite
    )
    if multiprocess:
        p = multiprocessing.Pool(n_cores)
        f_args = [(info_json_file, output_dir) for info_json_file in info_json_files]
        p.starmap(run_pe_w_params, f_args)
        p.close()
        p.join()
    else:
        for info_json_file in info_json_files:
            run_pe_w_params(info_json_file, output_dir)

    # output the parameters used
    params = {
        "gap_frac_cutoff": gap_frac_cutoff,
        "overwrite": overwrite,
        "json_directory": str(json_dir),
    }
    with open(output_dir / "params.json", "w") as f:
        json.dump(params, f, indent=4)


if __name__ == "__main__":
    main(
        json_dir=JSON_DIR,
        output_dir=OUTPUT_DIR,
        gap_frac_cutoff=GAP_FRAC_CUTOFF,
        overwrite=OVERWRITE,
        n_cores=N_CORES,
        multiprocess=True,
    )
