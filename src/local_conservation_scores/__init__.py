from .aln_asym_sum_of_pairs import main as score_aln_asym_sum_of_pairs
from .aln_property_entropy import main as score_aln_property_entropy


class conservation_score_methods:
    def __init__(self):
        self.aln_asym_sum_of_pairs = score_aln_asym_sum_of_pairs
        self.aln_property_entropy = score_aln_property_entropy
    
    def __getitem__(self, key):
        return getattr(self, key)
    
    def __setitem__(self, key, value):
        setattr(self, key, value)
    
    def __delitem__(self, key):
        delattr(self, key)

CONSERVATION_SCORE_METHODS = conservation_score_methods()

    
