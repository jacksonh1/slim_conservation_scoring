from pathlib import Path

import numpy as np
import pandas as pd
from Bio import Align, AlignIO, Seq, SeqIO
from pyprojroot import here

import local_seqtools.substitution_matrices as submats
import math

# %%
root = here()
data = root / "data"
disorder_matrix_folder = data / "disorder-matrix"
matrix_folder = data / "substitution_matrices"
matrix_folder.mkdir(exist_ok=True)
# %%

def prep_matrix_dfs(mat_df, mat_name, output_dir="./matrices"):
    mat_df = mat_df.copy()
    mat_df_new_diag = submats.matrixdf_diagonal_2_max_off_diagonal(mat_df)
    mat_df_norm = submats.normalize_matrix_df(mat_df)
    mat_df_new_diag_norm = submats.normalize_matrix_df(mat_df_new_diag)
    mat_df_sqrt_norm = submats.sqrt_normalize_matrix_df(mat_df)

    mat_df.to_csv(output_dir / f"{mat_name}.csv")
    mat_df_new_diag.to_csv(output_dir / f"{mat_name}_max_off_diagonal.csv")
    mat_df_norm.to_csv(output_dir / f"{mat_name}_norm.csv")
    mat_df_new_diag_norm.to_csv(output_dir / f"{mat_name}_max_off_diagonal_norm.csv")
    mat_df_sqrt_norm.to_csv(output_dir / f"{mat_name}_sqrt_norm.csv")


matrix_name = "BLOSUM62"
matBLOSUM62_df = submats.load_matrix_as_df(matrix_name)
prep_matrix_dfs(matBLOSUM62_df, matrix_name, output_dir=matrix_folder)

matrix_name = "EDSSMat50"
mat = Align.substitution_matrices.read(
    disorder_matrix_folder / "Matrices_and_Datasets/Matrices/EDSSMat50"
)
matEDSS50_df = submats.convert_matrix_array_2_df(mat)
prep_matrix_dfs(matEDSS50_df, matrix_name, output_dir=matrix_folder)
