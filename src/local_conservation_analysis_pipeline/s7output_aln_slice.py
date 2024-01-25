# %%

import json
from pathlib import Path

import alv
import local_conservation_analysis_pipeline.group_conservation_objects as group_tools
import local_env_variables.env_variables as env
from Bio import Align, AlignIO, SeqIO

json_file = "./conservation_analysis/2-9606_0_002f40/2-9606_0_002f40.json"
n_flanking_cols= 20

def save_colored_protein_msa_html(
    alignment: Align.MultipleSeqAlignment,
    output_html_file,
    color_scheme="clustal",
    max_id_length=None,
):
    """
    Save a colored view of a protein multiple sequence alignment (MSA) as HTML.

    Parameters:
    - alignment_file (str): Path to the input protein MSA file (e.g., in FASTA format).
    - output_html_file (str): Path to save the output HTML file with colored protein MSA.
    - color_scheme (str): Color scheme for the alignment. Options: "default", "clustal", "mismatch".

    Returns:
    None
    """
    if color_scheme == "clustal":
        with open(env.COLOR_MAP_FILES["clustal"], "r") as f:
            colors = json.load(f)
    elif color_scheme == "mismatch":
        colors = {
            "A": "blue",
            "R": "red",
            "N": "green",
            "D": "orange",
            "C": "purple",
            "Q": "cyan",
            "E": "yellow",
            "G": "brown",
            "H": "pink",
            "I": "gray",
            "L": "olive",
            "K": "darkgreen",
            "M": "darkblue",
            "F": "darkred",
            "P": "darkorange",
            "S": "darkmagenta",
            "T": "darkyellow",
            "W": "darkcyan",
            "Y": "darkgray",
            "V": "lightgray",
            "-": "lightgray",
        }
    else:
        # Default color scheme
        colors = {
            "A": "green",
            "R": "blue",
            "N": "purple",
            "D": "red",
            "C": "orange",
            "Q": "magenta",
            "E": "yellow",
            "G": "cyan",
            "H": "lightblue",
            "I": "brown",
            "L": "pink",
            "K": "gray",
            "M": "olive",
            "F": "darkgreen",
            "P": "darkblue",
            "S": "darkred",
            "T": "darkorange",
            "W": "darkmagenta",
            "Y": "darkyellow",
            "V": "darkcyan",
            "-": "lightgray",
        }

    html_content = "<html><head><style>pre {font-family: 'Courier New', monospace;}</style></head><body><pre>\n"
    max_id_length = max_id_length or max(len(record.id) for record in alignment)
    for record in alignment:
        formatted_id = record.id[:max_id_length].ljust(max_id_length)
        sequence_line = ""
        for symbol in record.seq:
            sequence_line += (
                f'<span style="color: {colors.get(symbol, "black")}">{symbol}</span>'
            )
        html_content += f"{formatted_id} {sequence_line}\n"
    html_content += "</pre></body></html>"
    with open(output_html_file, "w") as html_file:
        html_file.write(html_content)

def main(json_file, n_flanking_cols):
    og = group_tools.ConserGene(json_file)
    og.load_levels()
    for level in og.levels_passing_filters:
        lvl = og.level_objects[level]
        aln = lvl.aln
        start = max(0, lvl.hit_aln_start - n_flanking_cols)
        end = min(len(aln[0]), lvl.hit_aln_end + 1 + n_flanking_cols)
        aln_slice = slice(start, end)
        slice_file = Path(og.analysis_folder) / f"{og.reference_index}-{og.query_gene_id.replace(':','')}-{level}_aln_slice.html"
        og.add_item_to_lvl_orthogroup("aln_slice_file", str(slice_file), level)
        save_colored_protein_msa_html(aln[:, aln_slice], slice_file, color_scheme="clustal")
