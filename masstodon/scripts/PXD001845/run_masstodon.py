%load_ext autoreload
%autoreload 2
%load_ext line_profiler

from collections import Counter
from itertools  import islice
import json
import os
from os.path    import join as pjoin, exists as pexists, split as psplit
from time       import time


from masstodon.scripts.PXD001845.iter_folders import path
from masstodon.scripts.PXD001845.iter_folders import get_charge
from masstodon.scripts.PXD001845.iter_folders import iter_scans, non_modified_scans
from masstodon.scripts.PXD001845.hecks_custom_ptms import modify_fasta
from masstodon.masstodon import masstodon_single

class Suspicion(Exception):
    pass

with open("/Users/matteo/Projects/masstodon/data/PXD001845/csv2fasta.json", 'r') as f:
    csv2fasta = json.load(f)

def iter_data(path):
    for mz, intensity, charge, experiment in iter_scans(path):
        exp = experiment.split("/")[-2]
        if "AMB_pBora":
            yield mz, intensity, charge, experiment

def read_fasta(experiment, csv2fasta):
    folder = experiment.split('/')[-2]
    matching_csv = [f for f in csv2fasta if folder in f]
    if len(matching_csv) == 1:
        fasta = csv2fasta[matching_csv[0]]
        return fasta
    else:
        raise Suspicion("Suspicious to have more than one folder here.")


bad = []
if __name__ == "__main__":
    min_prob          = .8
    isotopic_coverage = .999
    std_cnt           = 3
    dump_path         = "/Users/matteo/Projects/masstodon/dumps/one_process"
    stop              = 100
    for mz, intensity, q, path in islice(iter_data(path), stop):
        exp = path.split('/')[-2:]
        print("Starting " + str(exp))
        try:
            fasta = read_fasta(path, csv2fasta)
            fasta, modifications = modify_fasta(fasta)
            M, timings = masstodon_single(mz, intensity, fasta, q, '',
                                          isotopic_coverage = isotopic_coverage,
                                          min_prob          = min_prob,
                                          std_cnt           = std_cnt,
                                          orbitrap          = True,
                                          get_timings       = True)
            local_dump_path = pjoin(dump_path, pjoin(*exp))
            if not pexists(local_dump_path):
                os.makedirs(local_dump_path)
            M.dump(local_dump_path, indent=4)
            M.write(local_dump_path)
            M.plotlygl(local_dump_path, show=False)
            with open(pjoin(local_dump_path, 'timings.json'), 'w') as h:
                json.dump(timings, h, indent=4)
            print("\tdone\n")
        except Exception as e:
            print(path)
            print(e)
            bad.append(dict(mz                  = mz,
                            intensity           = intensity,
                            fasta               = fasta,
                            q                   = q,
                            name                = '',
                            isotopic_coverage   = isotopic_coverage,
                            min_prob            = min_prob,
                            std_cnt             = std_cnt,
                            orbitrap            = True,
                            get_timings         = True))

# there is something wrong above:
previously_wrong = "20141202_AMB_Bora_10x_40MeOH_1FA_OT_120k_10uscans_795_ETciD_4ms_15SA_22precZ/1"


min_prob          = .8
isotopic_coverage = .999
std_cnt           = 3
dump_path         = "/Users/matteo/Projects/masstodon/dumps/one_process"
stop              = None
all_path = "/Users/matteo/Projects/masstodon/data/PXD001845/numpy_files/"

sol = []
# we match to both *1 and *10 ... idiot
for mz, intensity, q, path in islice(iter_data(all_path), stop):
    if previously_wrong in path:
        exp = path.split('/')[-2:]
        print("Starting " + str(exp))
        try:
            fasta = read_fasta(path, csv2fasta)
            fasta, modifications = modify_fasta(fasta)
            M, timings = masstodon_single(mz, intensity, fasta, q, '',
                                          isotopic_coverage = isotopic_coverage,
                                          min_prob          = min_prob,
                                          std_cnt           = std_cnt,
                                          orbitrap          = True,
                                          get_timings       = True)
            sol.append(M)
        except Exception as e:
            print(e)
    else:
        pass

M = sol[0]
M.plotlygl("dump")



visited_files = [path for mz, intensity, q, path in iter_data(all_path)]
visited_files[:40]