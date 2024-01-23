import json
from Bio import SeqIO
import local_env_variables.env_variables as env
import sys
import os
from pathlib import Path


def main(input_file, output_file):
    with open(input_file, 'r') as f:
        data = json.load(f)
    id_list = data['sequences_clustered_ldos']
    seqrecords = [env.ODB_DATABASE.data_all_seqrecords_dict[id] for id in id_list]
    aln_file = Path(output_file.replace('.fasta', '_clustered_ldos_aln.fasta'))
    data['alignment_clustered_ldos_file'] = str(aln_file.resolve())
    with open(input_file, 'w') as f:
        json.dump(data, f, indent=4)
    with open(output_file, 'w') as f:
        SeqIO.write(seqrecords, f, 'fasta')


if __name__ == '__main__':
    input_file = str(sys.argv[1])
    output_file = str(sys.argv[2])
    main(input_file, output_file)

