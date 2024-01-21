# %%
import copy
import json
import os
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import local_seqtools.general_utils as tools
import meta_tools
import numpy as np
from pyprojroot import here

import pandas as pd
from Bio import Align, AlignIO, Seq, SeqIO

# %load_ext autoreload
# %autoreload 2
# %%

OG_ANALYSIS_OUTPUT_FOLDER=meta_tools.OG_ANALYSIS_OUTPUT_FOLDER
HITS_FILE = meta_tools.HITS_FILE
import local_env_variables.env_variables_and_filepaths as fp
ROOT = fp.conservation_analysis_root

df = pd.read_csv(HITS_FILE.replace('.csv', '_original_reindexed.csv'))

with open('./s07_full_msa_annotation.json') as f:
    s07 = json.load(f)

# convert dict to dataframe
s07df = pd.DataFrame.from_dict(s07, orient='index')
s07df = s07df.reset_index(names=['reference_index'])
s07df['reference_index'] = s07df['reference_index'].astype(int)

dfm = pd.merge(df, s07df, on='reference_index', how='left')


meta_dict = copy.deepcopy(meta_tools.META_PARAMS)
failmap = {int(k): v for k, v in meta_dict['failed_id_dict'].items()}
note_map = {int(k): v for k, v in meta_dict['notes'].items()}


dfm['fail_reason'] = dfm['reference_index'].map(failmap)
dfm['notes'] = dfm['reference_index'].map(note_map)


# %%
ortholog_group_info_map = {}
for file in list(Path(OG_ANALYSIS_OUTPUT_FOLDER).glob('*/ortholog_group_info.json')):
    with open(file) as f:
        ortho_info = json.load(f)
    ref_index = ortho_info['reference_index']
    odbid = ortho_info['query_gene_id']
    odb_analysis_folder = Path(file).parent.resolve().relative_to(ROOT)
    if ortho_info['found hit']:
        hit_sequence = ortho_info['hit_sequence']
        hit_in_idr = ortho_info['hit_in_idr']
    else:
        hit_sequence = np.nan
        hit_in_idr = np.nan
    if 'multi_plot_filename' in ortho_info:
        multi_plot_filename = str(Path(ortho_info['multi_plot_filename']).name)
    else:
        multi_plot_filename = np.nan
    ortholog_group_info_map[ref_index] = [
        odbid, hit_sequence, hit_in_idr, multi_plot_filename, odb_analysis_folder
    ]



# %%


# map the dictionary to the dataframe
# make each element of the list a column
dfm[['orthoDB id (can search on orthoDB website)', 'hit sequence used for conservation', 'hit in idr (by custom iupred)', 'multi_plot_image_filename', 'odb_analysis_folder']] = dfm['reference_index'].map(ortholog_group_info_map).apply(pd.Series)

# %%

dfm.to_csv(HITS_FILE.replace('.csv','-ANALYZED.csv'), index=False)







# %%
