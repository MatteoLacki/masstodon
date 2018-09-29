from collections import Counter, defaultdict
import json
from os import listdir
from os.path import join as pjoin, split as psplit, exists as pexists
import pandas as pd
import numpy as np
from sys import platform

from masstodon.data.ptms import ptms

class WrongSystem(Exception):
    pass

def get_folder2ptms():
    if platform == "darwin":
        # check if you have your latte on skimmed soya milk with you
        datapath = "/Users/matteo/Projects/masstodon/data/PXD001845/csv_files"
    elif platform == "linux":
        # check if you have long dirty hair
        data_path = "/home/matteo/masstodon/review_answer/csv_files"
    elif "win" in platform:
        # don't check anything. no use.
        data_path = "C:/"
        raise WrongSystem(":D")
    else:
        raise WrongSystem("Path not specified correctly.")
    csvs = listdir(datapath)
    folder2ptms = defaultdict(list)
    for csv in csvs:
        csv = csv.split(".")[0]
        folder, mods = csv.split("XlinkX")
        folder = folder[:-7] 
        mods = mods.split("_")[1:]
        parsed_mods = {}
        for m in mods:
            aa = m[1]
            pos = int(m[2:]) - 1
            parsed_mods[pos+1] = {"C_alpha": ptms["phosphorylation"].copy()}
        if parsed_mods:
            folder2ptms[folder].append(parsed_mods)
    return folder2ptms