# %%
import pairk
import conservation_scores.tools.pairwise_tools as pw_tools


mats = pw_tools.import_pairwise_matrices(
    "./2-Vertebrata-frag_pairwise_gapless_embedding_lf5_rf5.json"
)

pairk_mats = pairk.PairkAln.from_file("./2-Vertebrata-pairk_aln_embedding_lf5_rf5.json")

# %%


# function that compares two dataframes
def compare_dfs(df1, df2):
    # print(df1.shape)
    # print(df2.shape)
    assert df1.shape == df2.shape
    assert df1.columns.tolist() == df2.columns.tolist()
    assert df1.index.tolist() == df2.index.tolist()
    assert (df1 == df2).all().all()
    print("Dataframes are equal")


def compare_files(old_method_file, new_method_file):
    mats = pw_tools.import_pairwise_matrices(old_method_file)
    pairk_mats = pairk.PairkAln.from_file(new_method_file)
    compare_dfs(
        mats["score_dataframe"].rename(columns={"reference_kmer": "query_kmer"}),
        pairk_mats.score_matrix,
    )
    compare_dfs(
        mats["subseq_dataframe"].rename(columns={"reference_kmer": "query_kmer"}),
        pairk_mats.orthokmer_matrix,
    )
    compare_dfs(
        mats["position_dataframe"].rename(columns={"reference_kmer": "query_kmer"}),
        pairk_mats.position_matrix,
    )


# %%
compare_dfs(
    mats["score_dataframe"].rename(columns={"reference_kmer": "query_kmer"}),
    pairk_mats.score_matrix,
)
compare_dfs(
    mats["subseq_dataframe"].rename(columns={"reference_kmer": "query_kmer"}),
    pairk_mats.orthokmer_matrix,
)
compare_dfs(
    mats["position_dataframe"].rename(columns={"reference_kmer": "query_kmer"}),
    pairk_mats.position_matrix,
)

# %%
compare_files(
    "./2-Vertebrata-frag_pairwise_gapless_embedding_lf5_rf5.json",
    "./2-Vertebrata-pairk_aln_embedding_lf5_rf5.json",
)
compare_files(
    "./3-Vertebrata-frag_pairwise_gapless_embedding_lf5_rf5.json",
    "./3-Vertebrata-pairk_aln_embedding_lf5_rf5.json",
)
compare_files(
    "./2-Vertebrata-fragpair_gapless_lf5_rf5_edssmat50.json",
    "./2-Vertebrata-pairk_aln_lf5_rf5_edssmat50.json",
)
compare_files(
    "./3-Vertebrata-fragpair_gapless_lf5_rf5_edssmat50.json",
    "./3-Vertebrata-pairk_aln_lf5_rf5_edssmat50.json",
)
compare_files(
    "./3-Vertebrata-fragpair_gapless_needleman_lf5_rf5_edssmat50.json",
    "./3-Vertebrata-pairk_aln_needleman_lf5_rf5_edssmat50.json",
)
compare_files(
    "./2-Vertebrata-fragpair_gapless_needleman_lf5_rf5_edssmat50.json",
    "./2-Vertebrata-pairk_aln_needleman_lf5_rf5_edssmat50.json",
)


# %%


# %%

# %%
pairk_mats.score_matrix
mats["subseq_dataframe"] == pairk_mats.score_matrix

mats.keys()
(
    mats["score_dataframe"].rename(columns={"reference_kmer": "query_kmer"})
    == pairk_mats.score_matrix
).all().all()
