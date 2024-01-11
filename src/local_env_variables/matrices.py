
import pandas as pd
import local_env_variables.env_variables as env

matrices = env.MATRIX_DIR


MATRIX_DF_DICT = {
    'BLOSUM62': matrices / "BLOSUM62.csv",
    'BLOSUM62_norm': matrices / "BLOSUM62_norm.csv",
    'BLOSUM62_max_off_diagonal': matrices / "BLOSUM62_max_off_diagonal.csv",
    'BLOSUM62_max_off_diagonal_norm': matrices / "BLOSUM62_max_off_diagonal_norm.csv",
    'EDSSMat50': matrices / "EDSSMat50.csv",
    'EDSSMat50_norm': matrices / "EDSSMat50_norm.csv",
    'EDSSMat50_row_norm': matrices / "EDSSMat50_row_norm.csv",
    'EDSSMat50_max_off_diagonal': matrices / "EDSSMat50_max_off_diagonal.csv",
    'EDSSMat50_max_off_diagonal_norm': matrices / "EDSSMat50_max_off_diagonal_norm.csv",
    'grantham_similarity_norm': matrices / "grantham_similarity_norm.csv",
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
    mat_df = pd.read_csv(MATRIX_DF_DICT[matrix_name], index_col=0)
    return mat_df

