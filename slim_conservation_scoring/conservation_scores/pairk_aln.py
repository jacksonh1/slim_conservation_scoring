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
    """run pairk alignment from a fasta alignment file

    Parameters
    ----------
    input_alignment_file : str | Path
        input fasta alignment file
    output_file : str | Path
        output file for pairk alignment results in json format
    reference_id : str
        the sequence id of the query sequence
    k : int
        the value of k for pairk alignment
    idr_aln_st : int
        start position of the IDR in the alignment
    idr_aln_end : int
        end position of the IDR in the alignment
    overwrite : bool, optional
        whether or not to overwrite the output file if it exists, by default False
    matrix_name : str, optional
        the name of the scoring matrix to use for pairk.pairk_alignment, by default "EDSSMat50". The available matrices can be viewed by importing pairk and running the function `pairk.print_available_matrices()`.
    """
    # check if output file exists
    if Path(output_file).exists() and not overwrite:
        print(f"{output_file} exists and overwrite is False")
        print("using file that is already there...")
        return
    idr_dict = pairk.utilities.fasta_MSA_to_idr_dict(
        input_alignment_file, idr_aln_st, idr_aln_end
    )
    alnres = pairk.pairk_alignment(
        idr_dict,
        query_id=reference_id,
        k=k,
        matrix_name=matrix_name,
    )
    alnres.write_to_file(output_file)
