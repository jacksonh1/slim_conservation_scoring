# %%
import local_seqtools.substitution_matrices as submats
import pandas as pd

df = pd.read_csv('./grantham.tsv', sep='\t', index_col=0)
df.to_csv('../grantham.csv')


df_norm = submats.normalize_matrix_df(df)
df_norm.to_csv('./grantham_norm.csv')
df_sim_norm = 1-df_norm
df_sim_norm.to_csv('./grantham_similarity_norm.csv')
# %%