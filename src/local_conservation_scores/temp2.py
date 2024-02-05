from pathlib import Path

import numpy as np
import pandas as pd
from Bio import Align, AlignIO, Seq, SeqIO

import local_conservation_analysis_pipeline.group_conservation_objects as group_tools
from local_conservation_score_tools import capra_singh_2007_scores, score_tools
from local_env_variables import matrices
from local_seqtools import alignment_tools as aln_tools
from local_seqtools import general_utils as tools
from local_seqtools import jch_alignment as jch_aln

# ==============================================================================
# //
# ==============================================================================
