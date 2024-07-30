# %%
from pathlib import Path
import pairk

# %%


def main(
    input_alignment_file: str | Path,
    output_file: str | Path,
    reference_id: str,
    k: int,
    idr_aln_st: int,
    idr_aln_end: int,
    overwrite: bool = False,
    matrix_name: str = "grantham_similarity_norm",
    **kwargs,
):
    # check if output file exists
    if Path(output_file).exists() and not overwrite:
        print(f"{output_file} exists and overwrite is False")
        print("using file that is already there...")
        return
    aligner = pairk.make_aligner(matrix_name)
    idr_dict = pairk.utilities.fasta_MSA_to_idr_dict(
        input_alignment_file, idr_aln_st, idr_aln_end
    )
    alnres = pairk.pairk_alignment_needleman(
        idr_dict,
        query_id=reference_id,
        k=k,
        aligner=aligner,
    )
    alnres.write_to_file(output_file)
