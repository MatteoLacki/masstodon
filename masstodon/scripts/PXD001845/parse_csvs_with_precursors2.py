    from collections import Counter, defaultdict
import json
from os import listdir
from os.path import join as pjoin, split as psplit, exists as pexists
import pandas as pd
import numpy as np
from sys import platform

from masstodon.data.ptms import ptms

# csv = "20141202_AMB_Bora_10x_40MeOH_1FA_OT_120k_10uscans_920_ETD_8ms_19precZ_DC_L2_XlinkX.csv"

def get_folder2ptms(datapath):
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
    return dict(folder2ptms)