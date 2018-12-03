%load_ext autoreload
%autoreload 2

import os
from os import listdir
from os.path import join as pjoin
from collections import Counter
import numpy as np
import json
import re

from masstodon.read.mzml import read_mzml
from masstodon.plot.spectrum import plot_spectrum
from masstodon.scripts.PXD001845.parse_csvs_with_precursors2 import get_folder2ptms

path    = "/Users/matteo/Projects/masstodon/data/PXD001845/mzmls/"
csvpath = "/Users/matteo/Projects/masstodon/data/PXD001845/csv_files"
files   = [x for x in listdir(path) if ".mzML" in x]

# independent, multicore execution of masstodons
ETD = list(filter(lambda f: "ETD" in f, files))

# parameters for masstodon
min_prob          = .8
isotopic_coverage = .999
std_cnt           = 3
stop              = None
fasta             = "GAASMMGDVKESKMQITPETPGRIPVLNPFESPSDYSNLHEQTLASPSVFKSTKLPTPGKFRWSIDQLAVINPVEIDPEDIHRQALYLSHSRIDKDVEDKRQKAIEEFFTKDVIVPSPWTDHEGKQLSQCHSSKCTNINSDSPVGKKLTIHSEKSD"
folder2ptms       = get_folder2ptms(csvpath)

# this has to be enhanced for data other than ETD
get_q = lambda f: int(re.sub('[^0-9]','', f.split("_")[-1]))
# for f in ETD:
#     try:
#         print(get_q(f))
#     except ValueError:
#         print(f)

def iter_data(data=ETD):
    for f in data:
        q     = get_q(f)
        mods  = folder2ptms.get(f.replace('.mzML',''), [])
        f     = pjoin(path, f)
        for s in read_mzml(f):
            scan_no   = int(s['id'].split("=")[1])
            mz        = s['m/z array']
            intensity = s['intensity array']
            yield mz, intensity, q, fasta, mods

for x in iter_data():
    
