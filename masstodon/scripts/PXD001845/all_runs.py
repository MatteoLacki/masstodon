%load_ext autoreload
%autoreload 2
%load_ext line_profiler

from collections import Counter
from itertools   import islice
import json
import os
from os.path import join as pjoin, exists as pexists, split as psplit
from sys import platform

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
elif platform == "linux":
    # check if you have long dirty hair
    data_path = "/home/matteo/masstodon/review_answer/numpy_files/"
    dump_path = "/mnt/disk/masstodon/dumps/many_processes/"
elif "win" in platform:
    # don't check anything. no use.
    data_path = "C:/"
    raise WrongSystem(":D")
else:
    raise WrongSystem("Path not specified correctly.")


# running fittings
min_prob          = .8
isotopic_coverage = .999
std_cnt           = 3
stop              = 100
# for mz, intensity, q, path in islice(iter_data(path), stop):

mz, intensity, q, path =next(iter_scans(data_path))

def single_run(mz, intensity, q, path):
    exp, scan = path.split('/')[-2:]
    scan = int(scan)
    print("start: {exp}\tscan: {scan}".format(exp=exp, scan=scan))
    row = {"q":q, "exp":exp, "scan":scan}
    try:
        fasta = read_fasta(path, csv2fasta)
        fasta, modifications = modify_fasta(fasta)
        M, timings = masstodon_single(mz, intensity, fasta, q, '',
                                      isotopic_coverage = isotopic_coverage,
                                      min_prob          = min_prob,
                                      std_cnt           = std_cnt,
                                      orbitrap          = True,
                                      get_timings       = True)
        local_dump_path = pjoin(dump_path, pjoin(exp))
        if not pexists(local_dump_path):
            os.makedirs(local_dump_path)
        M.dump(local_dump_path, indent=4)
        M.write(local_dump_path)
        M.plotlygl(local_dump_path, show=False)
        with open(pjoin(local_dump_path, 'timings.json'), 'w') as h:
            json.dump(timings, h, indent=4)

        # errors
        row.update(M.imperator.errors())

        # timings
        row.update({"t_"+str(n): T for n,T in timings})
        
        # intensities
        for s in ('ETDorHTR', 'ETnoD_PTR_fragments', 'ETnoD_precursor', 'PTR_precursor'):
            row["cz_simple."+str(s)] = int(M.cz_simple.intensities[s])
            row["cz."+str(s)]        = int(M.cz.intensities[s])

        row['success'] = True
    except Exception as e:
        row['success'] = False
    return row


# debug error calculation for a dp
model = sol.model
model.X
model.Y
model._coef

model.l1_abs()
model.res()

model.plot()
sol

cc = sol.cc
groups = M.imperator.groups




