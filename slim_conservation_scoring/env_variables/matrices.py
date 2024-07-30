
import os
from pathlib import Path

import dotenv
import pandas as pd
from Bio import Align

# src_dir = Path(os.path.dirname(__file__)).parent
# ROOT = src_dir.parent
# dotenv_path = ROOT / '.env'
# dotenv.load_dotenv(dotenv_path)
# MATRIX_DIR = ROOT / os.environ['SUBSTITUTION_MATRIX_DIR']

dotenv_path = os.path.join(os.path.dirname(__file__), '.env')
dotenv.load_dotenv(dotenv_path)
src_dir = Path(dotenv_path).parent.parent
MATRIX_DIR = src_dir / os.environ['SUBSTITUTION_MATRIX_DIR']



MATRIX_DF_DICT = {
    'BLOSUM62': MATRIX_DIR / "BLOSUM62.csv",
    'BLOSUM62_norm': MATRIX_DIR / "BLOSUM62_norm.csv",
    'BLOSUM62_max_off_diagonal': MATRIX_DIR / "BLOSUM62_max_off_diagonal.csv",
    'BLOSUM62_max_off_diagonal_norm': MATRIX_DIR / "BLOSUM62_max_off_diagonal_norm.csv",
    'EDSSMat50': MATRIX_DIR / "EDSSMat50.csv",
    'EDSSMat50_norm': MATRIX_DIR / "EDSSMat50_norm.csv",
    'EDSSMat50_row_norm': MATRIX_DIR / "EDSSMat50_row_norm.csv",
    'EDSSMat50_max_off_diagonal': MATRIX_DIR / "EDSSMat50_max_off_diagonal.csv",
    'EDSSMat50_max_off_diagonal_norm': MATRIX_DIR / "EDSSMat50_max_off_diagonal_norm.csv",
    'grantham_similarity_norm': MATRIX_DIR / "grantham_similarity_norm.csv",
}

def load_precomputed_matrix_df(matrix_name='BLOSUM62_max_off_diagonal_norm'):
    '''
    BLOSUM62\n
    BLOSUM62_norm\n
    BLOSUM62_row_norm\n
    BLOSUM62_max_off_diagonal\n
    BLOSUM62_max_off_diagonal_norm\n
    EDSSMat50\n
    EDSSMat50_norm\n
    EDSSMat50_row_norm\n
    EDSSMat50_max_off_diagonal\n
    EDSSMat50_max_off_diagonal_norm\n
    grantham_similarity_norm\n
    '''
    mat_df = pd.read_csv(MATRIX_DF_DICT[matrix_name], index_col=0)
    return mat_df

def matrix_df_to_dict(matrix_df: pd.DataFrame) -> dict[str, dict[str, float]]:
    '''
    It's much faster to access a dictionary than a pandas dataframe
    '''
    matrix_dict = {}
    for i, row in matrix_df.iterrows():
        if i in matrix_dict:
            raise ValueError(f"Duplicate key {i}")
        matrix_dict[i] = {}
        for j, val in row.items():
            matrix_dict[i][j] = val
    return matrix_dict

    
# these are compatible with biopython's Align aligner
ALIGNER_MATRIX_FILE_DICT = {
    'grantham_similarity': MATRIX_DIR / "grantham_similarity_normx100_aligner_compatible",
    'EDSSMat50': MATRIX_DIR / "EDSSMat50_aligner_compatible",
    'BLOSUM62': MATRIX_DIR / "BLOSUM62",
}

def load_matrix_for_aligner(matrix_name: str) -> Align.substitution_matrices.Array:
    if matrix_name in ALIGNER_MATRIX_FILE_DICT:
        mat = Align.substitution_matrices.read(ALIGNER_MATRIX_FILE_DICT[matrix_name])
    else:
        raise ValueError(f"Matrix {matrix_name} is not compatible with Aligner")
    return mat



# could do this in an object oriented way.
# The way that it is now, I basically import the library and treat it like an object
# so there's not much difference?