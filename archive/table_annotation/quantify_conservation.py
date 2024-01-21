import pandas as pd
import numpy as np
import local_env_variables.env_variables as env
from pathlib import Path
import yaml

with open('./params.yaml') as f:
    params 