import json
import os
import re
import sys
import traceback
from pathlib import Path

from Bio import Align, AlignIO, Seq, SeqIO
from pyprojroot import here

import local_env_variables.env_variables_and_filepaths as fp
import local_seqtools.general_utils as tools
from local_orthoDB_analysis_tools_v2 import \
    conservation_scoring_tools as score_tools

CONSERVATION_ROOT = fp.conservation_analysis_root


class level_info:
    """data about the hit sequence and the ortholog alignment for 1 phylogenetic level"""

    def __init__(self, aln_info_dict, root=CONSERVATION_ROOT, local_sequences=False, load_sequences=True):
        self.root = root
        self.info_dict = aln_info_dict
        self.level = self.info_dict["level"]
        self.UniprotID = self.info_dict["UniprotID"]
        self.reference_index = self.info_dict["reference_index"]
        self.name = self.info_dict["name"]
        self.sequence_files = self.info_dict["sequence_files"]
        # self.analysis_folder = root / self.info_dict["analysis_folder"]
        self.query_aligned_sequence_str = self.info_dict["query_aligned_sequence"]
        self.query_sequence_str = self.info_dict["query_sequence_str"]
        self.query_sequence_id_str = self.info_dict["query_sequence_id_str"]
        self.query_gene_id = self.info_dict["query_gene_id"]
        self.hit_sequence_str = self.info_dict["hit_sequence"]
        self.hit_start_position = self.info_dict["hit_start_position"]
        self.hit_end_position = self.info_dict["hit_end_position"]
        self.hit_alignment_start_position = self.info_dict[
            "hit_alignment_start_position"
        ]
        self.hit_alignment_end_position = self.info_dict["hit_alignment_end_position"]
        self.hit_in_idr = self.info_dict["hit_in_idr"]
        if self.hit_in_idr:
            self.idr_start = self.info_dict["idr_start"]
            self.idr_end = self.info_dict["idr_end"]
            self.idr_aln_start = self.info_dict["idr_aln_start"]
            self.idr_aln_end = self.info_dict["idr_aln_end"]
        if local_sequences:
            if 'local_sequence_files' not in self.info_dict:
                raise KeyError(f"`local_sequence_files` not in info_dict for {self.level} {self.name} {self.UniprotID} {self.query_sequence_id_str}")
            self.local_sequence_files = self.info_dict["local_sequence_files"]
            self.load_local_sequences()
        elif load_sequences:
            self.load_sequences()

    def __str__(self) -> str:
        s = f"""
        hit_score_class instance for {self.name} {self.UniprotID} {self.query_sequence_id_str}
        hit_in_idr: {self.hit_in_idr}
        level: {self.level}
        """
        return s

    def _load_alignment(self, file_path):
        self.loaded_alignment_file = file_path
        self.alignment_orthologs_ldo_clustered = AlignIO.read(file_path, "fasta")

    def load_sequences(self):
        self.sequences_full_OG = tools.import_fasta(
            self.root / self.sequence_files["fasta sequences - full OG"]
        )
        self.sequences_LDOs = tools.import_fasta(
            self.root / self.sequence_files["fasta sequences - OG LDOs"]
        )
        self.sequences_orthologs_ldo_clustered = tools.import_fasta(
            self.root / self.sequence_files["fasta sequences - OG LDOs cdhit"]
        )
        self.alignment_orthologs_ldo_clustered: AlignIO.MultipleSeqAlignment = (
            AlignIO.read(
                self.root / self.sequence_files["fasta alignment - OG LDO cdhit"],
                "fasta",
            )
        )

    def load_local_sequences(self):
        if "fasta sequences - full OG" in self.local_sequence_files:
            self.sequences_full_OG = tools.import_fasta(
                self.root / self.local_sequence_files["fasta sequences - full OG"]
            )
        if "fasta sequences - OG LDOs" in self.local_sequence_files:
            self.sequences_LDOs = tools.import_fasta(
                self.root / self.local_sequence_files["fasta sequences - OG LDOs"]
            )
        if "fasta sequences - OG LDOs cdhit" in self.local_sequence_files:
            self.sequences_orthologs_ldo_clustered = tools.import_fasta(
                self.root / self.local_sequence_files["fasta sequences - OG LDOs cdhit"]
            )
        if "fasta alignment - OG LDO cdhit" in self.local_sequence_files:
            self.alignment_orthologs_ldo_clustered: AlignIO.MultipleSeqAlignment = (
                AlignIO.read(
                    self.root
                    / self.local_sequence_files["fasta alignment - OG LDO cdhit"],
                    "fasta",
                )
            )

    def _slice_alignment(
        self,
        slice_start_position=None,
        slice_end_position=None,
    ):
        """
        Parameters
        ----------
        slice_start_position : int, optional
            start position of alignment slice. If not provided, the function uses `self.hit_alignment_start_position`. by default None
        slice_end_position : int, optional
            end position of alignment slice. If not provided, the function uses `self.hit_alignment_end_position+1`, by default None
        """
        if slice_start_position is None:
            slice_start_position = self.hit_alignment_start_position
        if slice_end_position is None:
            slice_end_position = self.hit_alignment_end_position + 1
        alignment_slice = self.alignment_orthologs_ldo_clustered[
            :, slice_start_position:slice_end_position
        ]
        query_seqrecord_in_alignment_slice = [
            i for i in alignment_slice if self.query_sequence_id_str in i.id
        ][0]
        # query_sequence_str_in_alignment_slice = str(self.query_seqrecord_in_alignment_slice.seq)
        return alignment_slice, query_seqrecord_in_alignment_slice


class level_score_class(level_info):
    def __init__(
        self, aln_info_dict, root=CONSERVATION_ROOT, local_sequences=False, score_key=None
    ):
        super().__init__(aln_info_dict, root=root, local_sequences=local_sequences)
        self.scores = None
        self.score_mask = None
        self.z_score_dict = None
        self.z_score_failure_reason = None
        if score_key is not None:
            if score_key in self.info_dict:
                self.load_scores(score_key)
            else:
                raise KeyError(f"score_key {score_key} not found in info_dict. No scores for this one.")

    def load_scores(self, score_key="asym_valday_score_json"):
        self.score_json = self.root / self.info_dict[score_key]
        with open(self.score_json, "r") as f:
            self.score_dict = json.load(f)
        self.scores = self.score_dict["scores"]
        self.score_mask = self.score_dict["score_mask"]

    def calculate_z_scores_bg_region(self, bg_region=None, num_bg_scores_cutoff=20):
        if bg_region is None:
            if self.hit_in_idr:
                bg_region = [self.idr_aln_start, self.idr_aln_end]
            else:
                self.z_score_failure_reason = (
                    "No bg region specified and hit is not in IDR"
                )
                return None
        try:
            self.z_score_dict = score_tools.calculate_z_score_bg_region(
                self.scores,
                self.score_mask,
                bg_region,
                num_bg_scores_cutoff=num_bg_scores_cutoff,
            )
        except ValueError as e:
            # traceback.print_exc()
            print(f"ValueError: {e}, {self.name} {self.UniprotID} {self.query_sequence_id_str}")
            self.z_score_failure_reason = str(e)
            return None
        if self.z_score_dict is not None:
            self.z_score_dict["bg_region"] = bg_region


class ortholog_group_info:
    """
    data about the entire ortholog group (all of the levels). Stored in 1 large json file for each protein
    loads and stores info about each level in the ortholog group as a level_info object
    """

    def __init__(self, ortholog_group_info_json, root=CONSERVATION_ROOT):
        self.root = root
        self.ortholog_group_info_json = Path(ortholog_group_info_json)
        with open(ortholog_group_info_json, "r") as f:
            self.info_dict = json.load(f)
        self.odb_folder = self.info_dict["odb_folder"]
        self.reference_index = self.info_dict["reference_index"]
        # self.analysis_folder = self.root / self.info_dict["analysis_folder"]
        self.UniprotID = self.info_dict["UniprotID"]
        self.name = self.info_dict["name"]
        self.query_gene_id = self.info_dict["query_gene_id"]
        self.levels = self.info_dict["levels"]
        self.found_hit = self.info_dict["found hit"]
        if self.found_hit:
            if "best orthogroup level" in self.info_dict:
                self.best_level = self.info_dict["best orthogroup level"]
            self.hit_sequence = self.info_dict["hit_sequence"]
            self.hit_in_idr = self.info_dict["hit_in_idr"]
            self.level_info_objects: dict[str, level_info] = {}
            self.level_score_objects: dict[str, level_score_class] = {}

    def __str__(self) -> str:
        s = f"""
        og_folder_class instance for {self.odb_folder}
        uniprot: {self.UniprotID}
        levels: {self.levels}
        file: {self.ortholog_group_info_json}"""
        return s

    def _overwrite_json(self):
        with open(self.ortholog_group_info_json, "w") as f:
            json.dump(self.info_dict, f, indent=4)

    def add_item_to_json(self, key, value, save_json=True):
        self.info_dict[key] = value
        if save_json:
            self._overwrite_json()

    def add_item_to_level_info(self, key, value, level, save_json=True):
        self.info_dict["level_info"][level][key] = value
        if save_json:
            self._overwrite_json()

    def load_level_info_object(self, level, local_sequences=False):
        lvl_info_obj = level_info(
            self.info_dict["level_info"][level],
            self.root,
            local_sequences=local_sequences,
        )
        return lvl_info_obj

    def load_all_level_info_objects(self, local_sequences=False):
        for level in self.levels:
            lvl_info_obj = self.load_level_info_object(
                level, local_sequences=local_sequences
            )
            self.level_info_objects[level] = lvl_info_obj

    def load_level_score_object_w_z_score(
        self,
        level,
        local_sequences=False,
        score_key="asym_valday_score_json",
        **z_score_params,
    ):
        try:
            og_scores_o = level_score_class(
                self.info_dict["level_info"][level],
                self.root,
                local_sequences=local_sequences,
                score_key=score_key,
            )
            og_scores_o.calculate_z_scores_bg_region(**z_score_params)
        except KeyError as e:
            # traceback.print_exc()
            # print(f"KeyError: {e}")
            raise e
            og_scores_o = None
        return og_scores_o

    def load_all_level_score_objects_with_z_score(
        self,
        local_sequences=False,
        score_key="asym_valday_score_json",
        **z_score_params,
    ):
        for level in self.levels:
            og_scores_o = self.load_level_score_object_w_z_score(
                level,
                local_sequences=local_sequences,
                score_key=score_key,
                **z_score_params,
            )
            self.level_score_objects[level] = og_scores_o


# class ortholog_group_info_scores(ortholog_group_info):

#     def __init__(self, ortholog_group_info_json, root=here()):
#         super().__init__(ortholog_group_info_json, root=root)
#         self.score_files = self.info_dict["score_files"]
#         if self.found_hit:
#             self.hit_score_objects: dict[str,hit_score_class.hit_score_class] = {}

#     def load_hit_score_object(self, level):
#         og_scores_o = hit_score_class.hit_score_class(self.root / self.score_files[level], self.root)
#         return og_scores_o

#     def load_all_hit_score_objects_with_z_score(self, **z_score_params):
#         for level in self.levels:
#             og_scores_o = self.load_hit_score_object(level)
#             og_scores_o.calculate_z_scores_bg_region(**z_score_params)
#             self.hit_score_objects[level] = og_scores_o


# class ortholog_group_info_scores_local_alignment(ortholog_group_info):

#     def __init__(self, ortholog_group_info_json, root=here()):
#         super().__init__(ortholog_group_info_json, root=root)
#         self.local_sequence_files = self.info_dict["local_sequence_files"]
#         if self.found_hit:
#             self.hit_score_objects: dict[str,hit_score_class.hit_score_class_local_alignment] = {}
#             self.aln_info_objects: dict[str,hit_score_class.local_aligment_info_from_score_json] = {}

#     def load_hit_score_object(self, level) -> hit_score_class.hit_score_class_local_alignment:
#         og_scores_o = hit_score_class.hit_score_class_local_alignment(self.root / self.score_files[level], self.root)
#         return og_scores_o

#     def load_all_hit_score_objects_with_z_score(self, **z_score_params):
#         for level in self.local_sequence_files.keys():
#             og_scores_o = self.load_hit_score_object(level)
#             og_scores_o.calculate_z_scores_bg_region(**z_score_params)
#             alignment_file = self.root / self.local_sequence_files[level]
#             og_scores_o.load_local_alignment(alignment_file)
#             self.hit_score_objects[level] = og_scores_o

#     def load_aln_info_object_from_score_file(self, level) -> hit_score_class.local_aligment_info_from_score_json:
#         aln_info_o = hit_score_class.local_aligment_info_from_score_json(self.root / self.score_files[level], self.root)
#         alignment_file = self.root / self.local_sequence_files[level]
#         aln_info_o.load_local_alignment(alignment_file)
#         return aln_info_o

#     def load_all_local_aln_info_objects_from_score_files(self):
#         for level in self.local_sequence_files.keys():
#             self.aln_info_objects[level] = self.load_aln_info_object_from_score_file(level)


# class level_info_with_sequences(level_info):

#     def __init__(self, aln_info_dict, root=here()):
#         super().__init__(aln_info_dict, root=root)
#         self.sequences_full_OG = tools.import_fasta(
#             root / self.sequence_files["fasta sequences - full OG"]
#         )
#         self.sequences_LDOs = tools.import_fasta(
#             root / self.sequence_files["fasta sequences - OG LDOs"]
#         )
#         self.sequences_orthologs_ldo_clustered = tools.import_fasta(
#             root / self.sequence_files["fasta sequences - OG LDOs cdhit"]
#         )
#         self.alignment_orthologs_ldo_clustered: AlignIO.MultipleSeqAlignment = AlignIO.read(
#             root / self.sequence_files["fasta alignment - OG LDO cdhit"],
#             "fasta",
#         )


# class level_info_with_sequences_local(level_info):

#     def __init__(self, aln_info_dict, root=here()):
#         super().__init__(aln_info_dict, root=root)
#         self.local_sequence_files = self.info_dict["local_sequence_files"]
#         if 'fasta sequences - full OG' in self.local_sequence_files:
#             self.sequences_full_OG = tools.import_fasta(
#                 root / self.local_sequence_files["fasta sequences - full OG"]
#             )
#         if 'fasta sequences - OG LDOs' in self.local_sequence_files:
#             self.sequences_LDOs = tools.import_fasta(
#                 root / self.local_sequence_files["fasta sequences - OG LDOs"]
#             )
#         if 'fasta sequences - OG LDOs cdhit' in self.local_sequence_files:
#             self.sequences_orthologs_ldo_clustered = tools.import_fasta(
#                 root / self.sequence_files["fasta sequences - OG LDOs cdhit"]
#             )
#         if 'fasta alignment - OG LDO cdhit' in self.local_sequence_files:
#             self.alignment_orthologs_ldo_clustered: AlignIO.MultipleSeqAlignment = AlignIO.read(
#                 root / self.sequence_files["fasta alignment - OG LDO cdhit"],
#                 "fasta",
#             )
