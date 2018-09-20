%load_ext autoreload
%autoreload 2
%load_ext line_profiler

from masstodon.read.txt import spectrum_from_txt
from masstodon.plot.spectrum import plot_spectrum

import numpy as np

datapath = "/Users/matteo/Projects/masstodon/data/Belgian/2013_05_01/SUBP/wave_height_1.5/FRL-010513-SUBP-WH 1,5-WV 10.txt"
mz, intensity = spectrum_from_txt(datapath)

# plot_spectrum(mz, intensity)

x = mz
w = intensity

threshold = 0.05

x - threshold
x + threshold

from math import inf

def smart_mz_iter(mz):
    MZ = iter(list(mz))
    prev_mz = -inf
    curr_mz = next(MZ)
    next_mz = next(MZ)
    yield prev_mz, curr_mz, next_mz
    for next_mz in MZ:
        prev_mz, curr_mz = curr_mz, next_mz
        yield prev_mz, curr_mz, next_mz
    prev_mz, curr_mz = curr_mz, next_mz
    yield prev_mz, curr_mz, inf

# sort!!!
l = []
r = []
for p, c, n in smart_mz_iter(mz):
    l.append(c-min((c-p)/2.0,threshold))
    r.append(c+min((n-p)/2.0,threshold))
l = np.array(l)
r = np.array(r)


SimpleGroups(l, r, x, w)