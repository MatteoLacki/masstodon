%load_ext autoreload
%autoreload 2
%load_ext line_profiler

from collections import Counter
import json
from os import listdir
from os.path import join as pjoin, split as psplit, exists as pexists
import pandas as pd
import numpy as np

datapath = "/Users/matteo/Projects/masstodon/data/PXD001845/csv_files"
fasta = "GAASMMGDVKESKMQITPETPGRIPVLNPFESPSDYSNLHEQTLASPSVFKSTKLPTPGKFRWSIDQLAVINPVEIDPEDIHRQALYLSHSRIDKDVEDKRQKAIEEFFTKDVIVPSPWTDHEGKQLSQCHSSKCTNINSDSPVGKKLTIHSEKSD"

csvs = listdir(datapath)
fastas = []
bad = []
for one_csv in csvs:
    try:
        df = pd.read_csv(pjoin(datapath, one_csv),
                     usecols=['pep_seq'])
        x = df['pep_seq'][~pd.isnull(df['pep_seq'])]
        longest_seq = max(x, key=len)
        fastas.append(longest_seq)
    except TypeError as te:
        print(te)
        print(one_csv)
        bad.append(one_csv)
    except ValueError as ve:
        try:
            df = pd.read_csv(pjoin(datapath, one_csv))
            x = df.iloc[1:,0]
            x = x[~pd.isnull(x)]
            longest_seq = max(x, key=len)
            fastas.append(longest_seq)
        except Exception as e:
            print(e)

# nic w tym śmiesznego...
all("SEKS" == f[-4:] for f in fastas)

# we add the bloody D at the end!
fastas = np.array(fastas)
fastas = np.core.defchararray.add(fastas, "D")

# there are only 9 types of different precursors here.
Counter(fastas)

# I have to make some sort of bloody dictionary.
csv2fasta = dict(zip(csvs, fastas))
with open("/Users/matteo/Projects/masstodon/data/PXD001845/csv2fasta.json", 'w') as f:
    json.dump(csv2fasta, f, indent=4)

# compare file names with csv2fasta keys...
file_names = listdir("/Users/matteo/Projects/masstodon/data/PXD001845/numpy_files/")

for fn in file_names:
    print(fn in csv2fasta)