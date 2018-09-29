import json
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

with open("/Users/matteo/Projects/masstodon/AC_plots/ETD_fit_stats.json", "r") as f:
    etd = json.load(f)

etd = pd.DataFrame(etd)
etd.to_csv('/Users/matteo/Projects/masstodon/AC_plots/ETD_fit_stats.csv')

def iter_estimates():
    for estimates, i in zip(etd.estimates, etd.index):
        for e in estimates:
            if e[0] != 'formula':
                out = {"i": i, 
                       "formula": e[0],
                       "q": e[1],
                       "g": e[2],
                       "intensity":e[4]}
                if len(e) == 6:
                    out["name"] = e[5]
                yield out
our_esitmates = pd.DataFrame(iter_estimates())
our_esitmates.to_csv("/Users/matteo/Projects/masstodon/AC_plots/estimates.csv", index=False)


# export map folder to csv
from collections import Counter, defaultdict
import json
from os import listdir
from os.path import join as pjoin, split as psplit, exists as pexists
import pandas as pd
import numpy as np
from masstodon.data.ptms import ptms


datapath = "/Users/matteo/Projects/masstodon/data/PXD001845/csv_files"
csvs = listdir(datapath)
folder2ptms = defaultdict(list)

def csv_folders():
    for csv in csvs:
        csv = csv.split(".")[0]
        folder, mods = csv.split("XlinkX")
        folder = folder[:-7] 
        yield csv, folder

pd.DataFrame(csv_folders()).to_csv("/Users/matteo/Projects/masstodon/AC_plots/csv2folders.csv", index=False)
# good



