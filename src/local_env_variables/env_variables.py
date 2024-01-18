import copy
import os
from pathlib import Path

import dotenv
import pandas as pd
# from attrs import define, field
from attrs import frozen
from Bio import SeqIO

dotenv.load_dotenv()

# get the root directory of the project (where the .env file is) and set it as the root
# This should be the root even if the script is run from a different directory
ROOT = Path(dotenv.find_dotenv()).parent

IUPRED_DIR = Path(os.environ['IUPRED2A_LIB_DIR'])

# use the root to convert to the correct absolute path
# MATRIX_DIR = Path(os.environ['SUBSTITUTION_MATRIX_DIR']).resolve().relative_to(ROOT).resolve()
MATRIX_DIR = ROOT / Path(os.environ['SUBSTITUTION_MATRIX_DIR'])
DATABASE_DIR = ROOT / os.environ['ORTHOGROUP_DATA_DIR']

phylogeny_lvl_ordering_file = ROOT / os.environ['PHYLOGENY_LVL_ORDERING']

with open(phylogeny_lvl_ordering_file, 'r') as f:
    PHYLOGENY_LVL_ORDERING = [l.strip() for l in f.readlines()]

COLOR_MAP_FILES = {
    'clustal': ROOT / 'data' / 'color_maps' / 'clustal_hex_map.json',
}
