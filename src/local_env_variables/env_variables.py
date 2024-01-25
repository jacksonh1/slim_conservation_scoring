import copy
import os
from pathlib import Path

import dotenv
import pandas as pd
# from attrs import define, field
from attrs import frozen
from Bio import SeqIO

# dotenv_path = os.path.join(os.path.dirname(__file__), '.env')
src_dir = Path(os.path.dirname(__file__)).parent
ROOT = src_dir.parent
dotenv_path = ROOT / '.env'
dotenv.load_dotenv(dotenv_path)

# this should be absolute path
IUPRED_DIR = os.environ['IUPRED2A_LIB_DIR']

phylogeny_lvl_ordering_file = ROOT / os.environ['PHYLOGENY_LVL_ORDERING']
with open(phylogeny_lvl_ordering_file, 'r') as f:
    PHYLOGENY_LVL_ORDERING = [l.strip() for l in f.readlines()]

colormap_dir = ROOT / os.environ['COLORMAP_DIR']
COLOR_MAP_FILES = {
    'clustal': colormap_dir / 'clustal_hex_map.json',
}



