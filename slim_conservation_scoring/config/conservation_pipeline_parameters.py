import copy
from pathlib import Path
from typing import Any, Literal, Union

from attrs import asdict, define, field, validators


@define
class HitSequenceParams:
    hit_sequence_search_method: Literal["search", "given_positions"] = field(
        default="search",
        validator=validators.in_(["search", "given_positions"]),  # type: ignore
    )  # type: ignore
    longest_common_subsequence: bool = field(default=False, converter=bool)
    lcs_min_length: int = field(default=20, converter=int)
    target_hit_length: int = field(default=0, converter=int)


@define
class IdrParams:
    find_idrs: bool = field(default=True, converter=bool)
    idr_map_file: str | Path | None = field(default=None)
    iupred_cutoff: float = field(
        default=0.4,
        converter=float,
        validator=validators.and_(validators.le(1), validators.ge(0)),
    )
    gap_merge_threshold: int = field(default=10, converter=int)
    idr_min_length: int = field(default=8, converter=int)

    def __attrs_post_init__(self):
        if self.find_idrs:
            assert (
                self.idr_map_file is None
            ), "idr_map_file must not be provided if find_idrs is True. choose one"
        else:
            assert (
                self.idr_map_file is not None
            ), "idr_map_file must be provided if find_idrs is False"


@define
class FilterParams:
    min_num_orthos: int = field(default=20, converter=int)


# @define
# class ScoreMethod:
#     """
#     lflank and rflank are only used for pairwise scores.
#     """

#     score_key: str = field(default="aln_property_entropy")
#     score_function_name: str = field(default="aln_property_entropy")
#     score_kwargs: dict[str, Any] = field(factory=dict)
#     level: str | None = field(default=None)
#     lflank: int = field(default=0, converter=int)
#     rflank: int = field(default=0, converter=int)


@define
class MSAScoreMethod:
    """
    lflank and rflank are only used for pairwise scores.
    """

    score_key: str = field(default="aln_property_entropy")
    function_name: str = field(default="aln_property_entropy")
    function_params: dict[str, Any] = field(factory=dict)
    level: str | None = field(default=None)


@define
class PairKAlnMethod:
    score_key: str = field(default="pairk_aln_needleman_lf0_rf0")
    function_name: str = field(default="pairk_aln_needleman")
    function_params: dict[str, Any] = field(factory=dict)
    level: str | None = field(default=None)
    lflank: int = field(default=4, converter=int)
    rflank: int = field(default=4, converter=int)


@define
class EsmParams:
    processes: int = field(default=4, converter=int)
    threads: int = field(default=1, converter=int)
    device: Literal["cuda", "cpu"] = field(default="cuda")
    model_name: str = field(default="esm2_t33_650M_UR50D")


@define
class PairKEmbeddingAlnMethod:
    """
    Score methods that use embeddings. You may want to define different parameters
    for running anything that uses embeddings because of the high memory usage and
    option to use GPU.
    """

    score_key: str = field(default="pairk_aln_embedding_lf0_rf0")
    function_name: str = field(default="pairk_aln_embedding")
    function_params: dict[str, Any] = field(factory=dict)
    level: str | None = field(default=None)
    lflank: int = field(default=0, converter=int)
    rflank: int = field(default=0, converter=int)


@define
class PairKmerConservationParams:
    kmer_conservation_function_name: str = "pairk_conservation"
    columnwise_score_function_name: str = "shannon_entropy"
    bg_cutoff: int = field(
        default=50,
        converter=int,
    )
    bg_kmer_cutoff: int = field(
        default=10,
        converter=int,
    )


@define
class MultiLevelPlotParams:
    score_key: str = field(default="aln_property_entropy")
    num_bg_scores_cutoff: int = field(default=20, converter=int)
    score_type: Literal["score", "z_score"] = field(
        default="z_score",
        validator=validators.in_(["score", "z_score"]),  # type: ignore
    )  # type: ignore
    strip_gaps: bool = field(default=False, converter=bool)


@define
class TableAnnotationConf:
    """ """

    score_key_for_table: str = field(default="aln_property_entropy")
    motif_regex: str | None = field(default=None)
    levels: list[str] = field(factory=list)
    annotations: list[str] = field(
        default=[
            "json_file",
            # "multi_level_plot",
            # "hit_start_position",
            # "regex",
            # "regex_match",
            # "regex_match_stpos_in_hit",
            # "conservation_string",
            # "aln_slice_file",
            # "hit_scores",
            # "hit_mean_score",
            # "hit_z_scores",
            # "hit_mean_zscore",
            # "best mean z-score over 5 residue window",
            # "best mean z-score over 10 residue window",
        ]
    )


@define
class AlnSliceConf:
    """
    parameters for creating the slice of the MSA surrounding the hit sequence

    Parameters
    ----------
    n_flanking_aas : int
        number of columns in the alignment to flank the hit sequence with
    """

    n_flanking_aas: int = field(default=20, converter=int)
    whole_idr: bool = field(default=False, converter=bool)


@define
class PipelineParameters:
    """
    Parameters for the pipeline
    """

    database_filekey: str | Path
    table_file: str | Path
    clear_files: bool = field(default=True, converter=bool)
    new_index: bool = field(default=True, converter=bool)
    steps_to_run: list = field(
        default=["s1", "s2", "s3", "s4", "s5", "s6", "s7", "s8", "s9"]
    )
    msa_score_methods: list[MSAScoreMethod] = field(factory=list)
    pairk_aln_methods: list[PairKAlnMethod] = field(factory=list)
    pairk_embedding_aln_methods: list[PairKEmbeddingAlnMethod] = field(factory=list)
    # pairwise_score_methods: list[PairwiseScoreMethod] = field(factory=list)
    output_folder: str | Path = field(default="conservation_analysis")
    hit_sequence_params: HitSequenceParams = field(default=HitSequenceParams())
    idr_params: IdrParams = field(default=IdrParams())
    filter_params: FilterParams = field(default=FilterParams())
    multilevel_plot_params: MultiLevelPlotParams = field(default=MultiLevelPlotParams())
    aln_slice_params: AlnSliceConf = field(default=AlnSliceConf())
    table_annotation_params: TableAnnotationConf = field(default=TableAnnotationConf())
    clean_analysis_files: bool = field(default=False, converter=bool)
    esm_params: EsmParams = field(default=EsmParams())
    pairk_conservation_params: PairKmerConservationParams = field(
        default=PairKmerConservationParams()
    )

    @classmethod
    def from_dict(cls, d: dict[str, Any]):
        d = copy.deepcopy(d)
        return cls(
            hit_sequence_params=HitSequenceParams(**d.pop("hit_sequence_params", {})),
            idr_params=IdrParams(**d.pop("idr_params", {})),
            filter_params=FilterParams(**d.pop("filter_params", {})),
            msa_score_methods=[
                MSAScoreMethod(score_key=k, **v)
                for k, v in d.pop("msa_score_methods", {}).items()
            ],
            pairk_aln_methods=[
                PairKAlnMethod(score_key=k, **v)
                for k, v in d.pop("pairk_aln_methods", {}).items()
            ],
            pairk_embedding_aln_methods=[
                PairKEmbeddingAlnMethod(score_key=k, **v)
                for k, v in d.pop("pairk_embedding_aln_methods", {}).items()
            ],
            esm_params=EsmParams(**d.pop("esm_params", {})),
            multilevel_plot_params=MultiLevelPlotParams(
                **d.pop("multilevel_plot_params", {})
            ),
            table_annotation_params=TableAnnotationConf(
                **d.pop("table_annotation_params", {})
            ),
            aln_slice_params=AlnSliceConf(**d.pop("aln_slice_params", {})),
            pairk_conservation_params=PairKmerConservationParams(
                **d.pop("pairk_conservation_params", {})
            ),
            **d,
        )

    def print_params(self):
        for k, v in asdict(self).items():
            if isinstance(v, dict):
                print(f"{k}:")
                for k2, v2 in v.items():
                    print(" ", f"{k2}:", v2)
                continue
            if isinstance(v, list):
                print(f"{k}:")
                for v2 in v:
                    print(" ", v2)
                continue
            print(f"{k}:", v)

    def check_for_dup_scorekeys(self):
        all_score_keys = []
        if len(self.msa_score_methods) > 0:
            for scoremethod in self.msa_score_methods:
                all_score_keys.append(scoremethod.score_key)
        if len(self.pairk_aln_methods) > 0:
            for pairk_method in self.pairk_aln_methods:
                all_score_keys.append(pairk_method.score_key)
        if len(self.pairk_embedding_aln_methods) > 0:
            for pairk_method in self.pairk_embedding_aln_methods:
                all_score_keys.append(pairk_method.score_key)
        seen = set()
        dupes = set()
        for score_key in all_score_keys:
            if score_key in seen:
                dupes.add(score_key)
            else:
                seen.add(score_key)
        if len(dupes) > 0:
            raise ValueError(
                f"there are duplicate score keys in the configuration file: {dupes}"
            )

    def get_score_key_dict(self):
        score_key_dict = {}
        if len(self.msa_score_methods) > 0:
            for scoremethod in self.msa_score_methods:
                score_key_dict[scoremethod.score_key] = scoremethod
        if len(self.pairk_aln_methods) > 0:
            for pairk_method in self.pairk_aln_methods:
                score_key_dict[pairk_method.score_key] = pairk_method
        if len(self.pairk_embedding_aln_methods) > 0:
            for pairk_method in self.pairk_embedding_aln_methods:
                score_key_dict[pairk_method.score_key] = pairk_method
        return score_key_dict

    def __attrs_post_init__(self):
        self.check_for_dup_scorekeys()

    #     if not Path(self.table_file).exists():
    #         raise FileNotFoundError(f"table_file {self.table_file} does not exist")
    #     if not Path(self.database_filekey).exists():
    #         raise FileNotFoundError(f"database_filekey {self.database_filekey} does not exist")


# print(test["new_score_methods"])

# import yaml

# config_file = "../../src/data_processing/benchmark_processing/p3_run_conservation_pipeline/params.yaml"
# with open(config_file, "r") as f:
#     config_dict = yaml.safe_load(f)
# test = {
#     "table_file": "../../benchmark/benchmark_v2/p1_table/benchmark_table_renamed.csv",
#     "database_filekey": "../../benchmark/benchmark_v2/p2_orthogroups/orthogroups/database_key.json",
# }
# config_dict.update(test)
# config = PipelineParameters.from_dict(config_dict)
# config = PipelineParameters(database_filekey="test", table_file="test")
# for k, v in asdict(config).items():
#     if isinstance(v, dict):
#         print(k)
#         for k2, v2 in v.items():
#             print("    ", k2, v2)
#         continue
#     if isinstance(v, list):
#         print(k)
#         for v2 in v:
#             print("    ", v2)
#         continue
#     print(k, v)
# print(config.score_methods)
# test = {
#     "output_folder": "./ortholog_analysis",
#     "database_filekey": "../../data/example_orthogroup_database_merged_version/human_odb_groups/database_key.json",
#     "table_file": "../../examples/table_annotation/table.csv",
#     # "hit_sequence_params": {
#         # "hit_sequence_search_method": "search"
#     # },
#     "idr_params": {
#         "find_idrs": True,
#         # idr_map_file: "./idr_map.json"
#         "iupred_cutoff": 0.4,
#         "gap_merge_threshold": 10,
#         "idr_min_length": 8,
#     },
#     "filter_params": {
#         "min_num_orthos": 20
#     },

#     "precalculated_aln_conservation_score_keys": [
#         "aln_property_entropy",
#     ],
#     "new_score_methods": {
#         "aln_property_entropy":{
#             "matrix_name": "EDSSMat50_max_off_diagonal_norm",
#             "gap_frac_cutoff": 0.2
#         },
#     },
#     "clear_files": False
# }
# test2 = {
#     "output_folder": "./ortholog_analysis",
#     "database_filekey": "../../data/example_orthogroup_database/human_odb_groups/database_key.json",
#     "table_file": "./table.csv",
#     "hit_sequence_params": {
#         "hit_sequence_search_method": "search"
#     },
#     "idr_params": {
#         "find_idrs": True,
#         # idr_map_file: "./idr_map.json"
#         "iupred_cutoff": 0.4,
#         "gap_merge_threshold": 10,
#         "idr_min_length": 8,
#     },
#     "filter_params": {
#         "min_num_orthos": 20
#     },
#     "new_score_methods": {
#         "aln_property_entropy":{
#             "matrix_name": "EDSSMat50_max_off_diagonal_norm",
#             "gap_frac_cutoff": 0.2
#         },
#     },
#     "clear_files": False
# }

# config = PipelineParameters.from_dict(test)
# for k, v in asdict(config).items():
#     print(k, v)

# config.precalculated_aln_conservation_score_keys.append("test")


# config2 = PipelineParameters.from_dict(test)
# for k, v in asdict(config2).items():
#     print(k, v)


# config = PipelineParameters.from_dict(test2)
# for k, v in asdict(config).items():
# print(k, v)
