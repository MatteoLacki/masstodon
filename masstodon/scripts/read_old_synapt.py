%load_ext autoreload
%autoreload 2
%load_ext line_profiler

from masstodon.read.txt import spectrum_from_txt
from masstodon.plot.spectrum import plot_spectrum
from masstodon.spectrum.threshold import threshold_spec

import numpy as np

datapath = "/Users/matteo/Projects/masstodon/data/Belgian/2013_05_01/SUBP/wave_height_1.5/FRL-010513-SUBP-WH 1,5-WV 10.txt"
mz, intensity = spectrum_from_txt(datapath)
threshold = 0.05


ts = threshold_spec(mz, intensity, threshold)
ts.get_smallest_diff_digits()
