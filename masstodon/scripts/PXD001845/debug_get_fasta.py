# %load_ext autoreload
# %autoreload 2
# %load_ext line_profiler

from itertools import islice
import json
import os
from os.path import join as pjoin, exists as pexists, split as psplit
from sys import platform
import multiprocessing as mp
from time import time

# Project specific
from masstodon.scripts.PXD001845.iter_folders import get_charge, iter_scans, non_modified_scans
from masstodon.scripts.PXD001845.get_fasta import read_fasta
from masstodon.scripts.PXD001845.hecks_custom_ptms import modify_fasta
from masstodon.scripts.PXD001845.csv2fasta import csv2fasta

from masstodon.masstodon import masstodon_single


# data paths
class WrongSystem(Exception):
    pass

if platform == "darwin":
    # check if you have your latte on skimmed soya milk with you
    data_path = "/Users/matteo/Projects/masstodon/data/PXD001845/numpy_files/"
    dump_path = "/Users/matteo/Projects/masstodon/dumps/many_processes/"
    processes_no = 1


# running fittings
min_prob          = .8
isotopic_coverage = .999
std_cnt           = 3
stop              = 25

folder = '20141202_AMB_pBora_AurA_10x_40MeOH_1FA_OT_120k_10uscans_928_ETD_8ms_19precZ' 
matching_csv = [f for f in csv2fasta if folder in f]

matching_csv[0]
csv2fasta[matching_csv[0]]

modify_fasta()


scans = list(iter_scans(data_path))

# testing
"20141202_AMB_pBora_AurA_10x_40MeOH_1FA_OT_120k_10uscans_928_ETD_8ms_19precZ_DC_L2_XlinkX_PS64_pT149"
f = "GAASMMGDVKESKMQITPETPGRIPVLNPFESPSDYSNLHEQTLASPSVFKSTKLPTPGKFRWJIDQLAVINPVEIDPEDIHRQALYLSHSRIDKDVEDKRQKAIEEFFTKDVIVPSPWTDHEGKQLSQCHSSKCTNINSDSPVGKKLBIHSEKSD"

modify_fasta(f)
f[149]