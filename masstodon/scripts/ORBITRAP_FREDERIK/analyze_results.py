import json
import numpy as np
import pandas as pd
from collections import Counter

import matplotlib.pyplot as plt


with open("/Users/matteo/Projects/masstodon/data/Belgian/2015_07/UBI_ORBI/12plusRes.json", "r") as f:
    etd = json.load(f)

def names_the_same(etd):
    """Are all names of fragments the same?"""
    for d in etd:
        for i, e in enumerate(d["estimates"]):
            if i>0:
                names = e[5]
                yield all(x == names[0] for x in names)
assert all(names_the_same(etd))

def get_intensities():
    intensities = Counter()
    for d in etd:
        try:
            for i, e in enumerate(d["estimates"]):
                if i > 0: # 0th row contains csv headers.
                    intensity = e[4]
                    name = e[5][0]
                    intensities[(d["exp"], name)] += intensity
        except IndexError:
            print(d["estimates"])
    return intensities
intensities = pd.DataFrame((k[0],k[1],v) for k,v in get_intensities().items())
intensities.columns = ['exp', 'frag', 'intensity']
intensities.to_csv('/Users/matteo/Projects/masstodon/AC_plots/masstodon/data/december/fred_ubi_intensities.csv', index=False)

etd = pd.DataFrame(etd)
etd.to_csv('/Users/matteo/Projects/masstodon/AC_plots/masstodon/data/december/fred_ubi_ETD_fit_stats.csv')

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
                    out["name"] = e[5][0]
                yield out
our_esitmates = pd.DataFrame(iter_estimates())
our_esitmates.to_csv("/Users/matteo/Projects/masstodon/AC_plots/masstodon/data/december/fred_ubi_estimates.csv", index=False)
