from .aln_asym_sum_of_pairs import main as score_aln_asym_sum_of_pairs
from .aln_property_entropy import main as score_aln_property_entropy
from .aln_shannon_entropy import main as score_aln_shannon_entropy
from .pairk_aln_embedding import main as pairk_embedding
from .pairk_aln import main as _pairk_aln
from .pairk_aln_needleman import main as _pairk_aln_needleman
from functools import partial
from slim_conservation_scoring.conservation_scores.tools import capra_singh_2007_scores
from .pairk_conservation import pairk_conservation_from_json


class MSAScoreMethods:
    """
    MSA-based conservation score methods
    """

    def __init__(self):
        self.aln_asym_sum_of_pairs = score_aln_asym_sum_of_pairs
        self.aln_property_entropy = score_aln_property_entropy
        self.aln_shannon_entropy = score_aln_shannon_entropy
        # self.fragment_pairwise_gapless = score_fragment_pairwise_gapless

    def __getitem__(self, key):
        return getattr(self, key)

    def __setitem__(self, key, value):
        setattr(self, key, value)

    def __delitem__(self, key):
        delattr(self, key)


# CONSERVATION_SCORE_METHODS = conservation_score_methods()


class PairKmerAlnMethods:
    """
    Pairk alignment methods
    """

    def __init__(self):
        self.pairk_aln = _pairk_aln
        self.pairk_aln_needleman = _pairk_aln_needleman
        self.pairk_aln_embedding = pairk_embedding

    def __getitem__(self, key):
        return getattr(self, key)

    def __setitem__(self, key, value):
        setattr(self, key, value)

    def __delitem__(self, key):
        delattr(self, key)


class PairKmerConservationMethods:
    """
    Pairk conservation methods
    """

    def __init__(self):
        # self.pairwise_matrix_to_kmer_scores = matrix_json_2_pairwise_scores
        self.pairk_conservation = pairk_conservation_from_json

    def __getitem__(self, key):
        return getattr(self, key)

    def __setitem__(self, key, value):
        setattr(self, key, value)

    def __delitem__(self, key):
        delattr(self, key)


class ColumnwiseScoreMethods:
    """
    Columnwise conservation score methods
    """

    def __init__(self):
        self.shannon_entropy = capra_singh_2007_scores.shannon_entropy
        self.property_entropy = capra_singh_2007_scores.property_entropy

    def __getitem__(self, key):
        return getattr(self, key)

    def __setitem__(self, key, value):
        setattr(self, key, value)

    def __delitem__(self, key):
        delattr(self, key)
