import multiprocessing
import json
from Bio import SeqIO
import sys
import os
from pathlib import Path
import local_env_variables.env_variables as env

database_dir = Path('../../../data/example_orthogroup_database_snakemake/human_odb_groups/')
JSON_DIR = database_dir / 'info_jsons'
OUTPUT_DIR = database_dir / 'clustered_ldo_fastas'
N_CORES = multiprocessing.cpu_count()





def save_fasta_from_json(
    json_file,
    output_dir,
):
    with open(json_file, 'r') as f:
        data = json.load(f)
    output_filename = Path(json_file).name.replace('.json', '.fasta')
    output_dir = Path(output_dir)
    output_file = output_dir / output_filename
    id_list = data['sequences_clustered_ldos']
    seqrecords = [env.ODB_DATABASE.data_all_seqrecords_dict[i] for i in id_list]
    aln_file_name = output_filename.replace('.fasta', '_clustered_ldos_aln.fasta')
    aln_file = output_dir.parent / 'alignments' / aln_file_name
    data['alignment_clustered_ldos_file'] = str(aln_file.resolve())
    with open(json_file, 'w') as f:
        json.dump(data, f, indent=4)
    with open(output_file, 'w') as f:
        SeqIO.write(seqrecords, f, 'fasta')

def main(
    json_dir=JSON_DIR,
    output_dir=OUTPUT_DIR,
    n_cores=N_CORES,
):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    info_json_files = list(json_dir.glob('*.json'))
    p = multiprocessing.Pool(n_cores)
    f_args = [(json_file, output_dir) for json_file in info_json_files]
    p.starmap(save_fasta_from_json, f_args)
    p.close()
    p.join()


if __name__ == '__main__':
    main(
        json_dir=JSON_DIR,
        output_dir=OUTPUT_DIR,
        n_cores=N_CORES,
    )


