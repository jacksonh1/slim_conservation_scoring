import copy
import os
from pathlib import Path

import dotenv
import pandas as pd
# from attrs import define, field
from attrs import frozen
from Bio import SeqIO

dotenv_path = os.path.join(os.path.dirname(__file__), '.env')
dotenv.load_dotenv(dotenv_path)
src_dir = Path(dotenv_path).parent.parent

# this should be absolute path
IUPRED_DIR = os.environ['IUPRED2A_LIB_DIR']

phylogeny_lvl_ordering_file = src_dir / os.environ['PHYLOGENY_LVL_ORDERING']
with open(phylogeny_lvl_ordering_file, 'r') as f:
    PHYLOGENY_LVL_ORDERING = [l.strip() for l in f.readlines()]

colormap_dir = src_dir / os.environ['COLORMAP_DIR']
COLOR_MAP_FILES = {
    'clustal': colormap_dir / 'clustal_hex_map.json',
}

