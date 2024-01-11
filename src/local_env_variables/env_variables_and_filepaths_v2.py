from dataclasses import dataclass

import pandas as pd
from attrs import define, field
from Bio import Seq, SeqIO
from pyprojroot import here

# ==============================================================================
# // general env variables and filepaths
# ==============================================================================

root = here()
data = root / "data"
orthodb_query_groups_root = data / 'orthoDB_ortholog_group_runs'
conservation_analysis_root = root / "conservation_analysis_runs"

# ==============================================================================
# // important filepaths
# ==============================================================================

orthodb11v0 = data / "orthoDB"
all_human_seq_fasta = orthodb11v0 / "odb_proteome" / "all_human_proteins_in_odb.fasta"
all_human_seq_clustered_fasta = orthodb11v0 / "odb_proteome" / "all_human_proteins_in_odb_clustered_c1.fasta"
pipeline_script_example_folder = root / "pipeline_scripts" / "orthoDB_table_analysis"


# ==============================================================================
# // ODB filepaths
# ==============================================================================
@define
class orthoDB_files_object:
    # gene_xrefs_full_tsv: str = str(orthodb11v0 / "odb11v0_gene_xrefs.tab")
    # gene_xrefs_uniprot_tsv: str = str(orthodb11v0 / "odb11v0_gene_xrefs_grepd_Uniprot_JCH.tab")
    # gene_xrefs_uniprot_human_tsv: str = str(orthodb11v0 / "odb11v0_gene_xrefs_grepd_Uniprot_9606_JCH.tab")
    # gene_refs_full_tsv: str = str(orthodb11v0 / "odb11v0_genes.tab")
    # ortholog_groups_tsv: str = str(orthodb11v0 / "odb11v0_OGs.tab")
    all_seqs_fasta: str = str(orthodb11v0 / "odb11v0_all_og.fasta")
    all_seqs_sqlite: str = str(orthodb11v0 / "odb11v0_all_og.sqlite")
    gene_xrefs_full_sqlite = str(orthodb11v0 / "odb11v0_gene_xrefs.sqlite")
    gene_refs_full_sqlite: str = str(orthodb11v0 / "odb11v0_genes.sqlite")
    ogs_sqlite: str = str(orthodb11v0 / "odb11v0_OGs.sqlite")
    OG2genes_sqlite: str = str(orthodb11v0 / "odb11v0_OG2genes.sqlite")
    levels_tsv: str = str(orthodb11v0 / "odb11v0_levels.tab")
    species_tsv: str = str(orthodb11v0 / "odb11v0_species.tab")
    readme: str = str(orthodb11v0 / "README.txt")
    JCH_readme: str = str(orthodb11v0 / "readme_JCH.txt")


orthoDB_files = orthoDB_files_object()





# ==============================================================================
# // ODB filepaths
# ==============================================================================
DATABASE_FILES = {
    "all sequences - fasta": orthodb11v0 / "odb11v0_all_og.fasta",
    "all sequences - sqlite": orthodb11v0 / "odb11v0_all_og.sqlite",
    "gene xrefs full - tsv": orthodb11v0 / "odb11v0_gene_xrefs.tab",
    "gene xrefs uniprot - tsv": (
        orthodb11v0 / "odb11v0_gene_xrefs_grepd_Uniprot_JCH.tab"
    ),
    "gene xrefs uniprot human - tsv": (
        orthodb11v0 / "odb11v0_gene_xrefs_grepd_Uniprot_9606_JCH.tab"
    ),
    "gene refs full - tsv": orthodb11v0 / "odb11v0_genes.tab",
    "gene refs full - sqlite": orthodb11v0 / "odb11v0_genes.sqlite",
    "ortholog groups (OGs) - tsv": orthodb11v0 / "odb11v0_OGs.tab",
    "levels": orthodb11v0 / "odb11v0_levels.tab",
    "species": orthodb11v0 / "odb11v0_species.tab",
    "readme": orthodb11v0 / "README.txt",
    "JCH readme": orthodb11v0 / "readme_JCH.txt",
    "OG2genes - sqlite": orthodb11v0 / "odb11v0_OG2genes.sqlite",
    "OGs per human gene - json": (
        orthodb11v0 / "s01_OG_groups_for_each_human_gene.json"
    ),
    "gene refs human - tsv": (orthodb11v0 / "odb11v0_genes_grepd_9606_0_JCH.tab"),
}


DATABASE_FILE_DESCRIPTIONS = {
    "all sequences - fasta": "all of the sequences in the database, in fasta format. Id's are the orthoDB gene ids",
    "all sequences - sqlite": "sqlite database for all of the sequences in the database",
    "gene xrefs full - tsv": "full list of alternate ids for each orthoDB gene id. For example UniProt ids, etc.. They are additional mappings not in the main gene refs file (like alternate uniprot ids). The full file is huge (19 GB)",
    "gene xrefs uniprot - tsv": "gene xrefs file, but filtered to include only UniProt ids. 2 GB",
    "gene xrefs uniprot human - tsv": "gene xrefs file, but filtered to include only UniProt ids for human genes. ~ 1 MB",
    "gene refs full - tsv": "full list of gene ids, their orthoDB species ids, other ids (like Uniprot), and a description. Very large file (8.1 GB)",
    "gene refs full - sqlite": "sqlite database for the gene refs full file",
    "ortholog groups (OGs) - tsv": "515 MB file. Has the orthoDB OG ID and the taxon level that the group is constructed at. Also has a name for the OG, that is named after the most common gene in the group.",
    "levels": "defines the OG taxon levels. Will correspond to the level in the OGs file. Also has the name of the level, and stats like the number of species below that level",
    "species": "has NCBI taxon id, orthoDB species id, and the name of the species",
    "readme": orthodb11v0 / "README.txt",
    "JCH readme": orthodb11v0 / "readme_JCH.txt",
    "OG2genes - sqlite": "sqlite database containing the different gene_ids for each OG_id"
}

# ==============================================================================
# // functions to print info to screen
# ==============================================================================
def print_database_filepaths():
    for k, v in DATABASE_FILES.items():
        print(f"- {k}\n    {v}")

def print_database_file_descriptions():
    for k, v in DATABASE_FILE_DESCRIPTIONS.items():
        print(f"- {k}\n    {v}")

def print_orthoDB_readme():
    with open(DATABASE_FILES["readme"], "r") as f:
        print(f.read())

# ==============================================================================
# // functions to load data
# ==============================================================================

def load_data_levels_df():
    data_levels_df = pd.read_csv(
        DATABASE_FILES["levels"],
        sep="\t",
        header=None,
        names=[
            "level NCBI tax id",
            "level name",
            "total non-redundant count of genes in all underneath clustered species",
            "total count of OGs built on it",
            "total non-redundant count of species underneath",
        ],
    )
    return data_levels_df


def load_data_uniprot_gene_key_dfs():
    """
    Relevant files:
        gene xrefs uniprot human - tsv
        gene refs human - tsv
        gene refs full - tsv
    """
    data_human_gene_key_df = pd.read_csv(
        DATABASE_FILES["gene refs human - tsv"],
        sep="\t",
        header=None,
        names=[
            "gene ID",
            "species ID",
            "source ID",
            "synonyms",
            "UniprotID",
            "Ensemble",
            "NCBI id/name",
            "description",
        ],
    )
    data_human_xref_gene_key_df = pd.read_csv(
        DATABASE_FILES["gene xrefs uniprot human - tsv"],
        sep="\t",
        header=None,
        names=["gene ID", "UniprotID", "DB name"],
    ).drop("DB name", axis=1)
    # return {"gene_key": data_human_gene_key_df, "xref_gene_key": data_human_xref_gene_key_df}
    return data_human_gene_key_df, data_human_xref_gene_key_df


def load_data_all_seqs():
    data_all_seqrecords_dict = SeqIO.index_db(
        str(DATABASE_FILES["all sequences - sqlite"]),
        str(DATABASE_FILES["all sequences - fasta"]),
        "fasta",
    )
    return data_all_seqrecords_dict


def load_data_og_df():
    data_og_df = pd.read_csv(
        DATABASE_FILES["ortholog groups (OGs) - tsv"],
        sep="\t",
        header=None,
        names=["OG id", "level NCBI tax id", "OG name"],
    )
    return data_og_df


def load_species_df():
    odb_species_df = pd.read_csv(
        DATABASE_FILES["species"],
        sep="\t",
        header=None,
        names=[
            "NCBI id",
            "species ID",
            "species name",
            "assembly ID",
            "n clustered genes",
            "n OGs",
            "mapping type",
        ],
    )
    return odb_species_df

def load_xref_uniprot():
    data_xref_uniprot_gene_key_df = pd.read_csv(
        DATABASE_FILES["gene xrefs uniprot - tsv"],
        sep="\t",
        header=None,
        names=["gene ID", "UniprotID", "DB name"],
    ).drop("DB name", axis=1)
    return data_xref_uniprot_gene_key_df

