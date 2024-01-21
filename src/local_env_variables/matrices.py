
import os
from pathlib import Path

import dotenv
import pandas as pd

dotenv.load_dotenv()
ROOT = Path(dotenv.find_dotenv()).parent
MATRIX_DIR = ROOT / Path(os.environ['SUBSTITUTION_MATRIX_DIR'])

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

def load_precomputed_matrix_df(matrix_name=None):
    '''
    BLOSUM62
    BLOSUM62_norm
    BLOSUM62_row_norm
    BLOSUM62_max_off_diagonal
    BLOSUM62_max_off_diagonal_norm
    EDSSMat50
    EDSSMat50_norm
    EDSSMat50_row_norm
    EDSSMat50_max_off_diagonal
    EDSSMat50_max_off_diagonal_norm
    grantham_similarity_norm
    '''
    if matrix_name is None:
        print("Available matrices:")
        for k in MATRIX_DF_DICT.keys():
            print(k)
        return
    mat_df = pd.read_csv(MATRIX_DF_DICT[matrix_name], index_col=0)
    return mat_df

