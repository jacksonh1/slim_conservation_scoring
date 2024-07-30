# %%
from pathlib import Path
import pairk


def main(
    input_alignment_file: str | Path,
    output_file: str | Path,
    reference_id: str,
    k: int,
    idr_aln_st: int,
    idr_aln_end: int,
    mod: pairk.ESM_Model,
    overwrite: bool = False,
    **kwargs,
):
    # check if output file exists
    if Path(output_file).exists() and not overwrite:
        print(f"{output_file} exists and overwrite is False")
        print("using file that is already there...")
        return
    idr_pos_map = pairk.utilities.fasta_MSA_to_idr_map(
        input_alignment_file, idr_aln_st, idr_aln_end
    )
    fl_seqdict_unaligned = pairk.utilities.fasta_MSA_to_unaligned_sequences(
        input_alignment_file
    )
    alnres = pairk.pairk_alignment_embedding_distance(
        full_length_sequence_dict=fl_seqdict_unaligned,
        idr_position_map=idr_pos_map,
        query_id=reference_id,
        k=k,
        mod=mod,
        **kwargs,
    )
    alnres.write_to_file(output_file)


# %%
# ==============================================================================
# // testing
# ==============================================================================
# import pipeline.group_conservation_objects as group_tools

# json_file = "../../benchmark/benchmark_multi_species/p3_conservation_analysis_pipeline/conservation_analysis/2-9606_0_004caa/2-9606_0_004caa.json"
# og = group_tools.ConserGene(json_file)
# lvlo=og.get_level_obj('Vertebrata')

# def pad_hit(seq: str, st_pos: int, end_pos: int, l_flank: int = 0, r_flank: int = 0):
#     st = max(0, st_pos - l_flank)
#     end = min(len(seq)-1, end_pos + r_flank)
#     return st, end, seq[st : end + 1]

# _,_,flanked_hit = pad_hit(og.query_sequence, og.hit_start_position, og.hit_end_position, 5, 5)
# k = len(flanked_hit)
# og.hit_sequence
# main(
#     input_alignment_file=lvlo.alignment_file,
#     reference_id=og.query_gene_id,
#     idr_aln_st=lvlo.idr_aln_start,
#     idr_aln_end=lvlo.idr_aln_end,
#     k=k,
#     output_file='test.json',
#     model_name='esm2_t33_650M_UR50D',
#     overwrite=True,
# )
