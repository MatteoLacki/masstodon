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

from masstodon.scripts.PXD001845.iter_folders import iter_scans
from masstodon.scripts.PXD001845.parse_csvs_with_precursors2 import get_folder2ptms
from masstodon.masstodon import masstodon_single, masstodon_batch


# data paths
class WrongSystem(Exception):
    pass

if platform == "darwin":
    # check if you have your latte on skimmed soya milk with you
    data_path = "/Users/matteo/Projects/masstodon/data/PXD001845/numpy_files/"
    dump_path = "/Users/matteo/Projects/masstodon/dumps/many_processes/"
    csvpath   = "/Users/matteo/Projects/masstodon/data/PXD001845/csv_files"
    processes_no = 1
elif platform == "linux":
    # check if you have long dirty hair
    data_path = "/home/matteo/masstodon/review_answer/numpy_files/"
    # dump_path = "/mnt/disk/masstodon/dumps/many_processes/"
    dump_path = "/mnt/disk/masstodon/dumps/no_charge_limits/"
    csvpath = "/home/matteo/masstodon/review_answer/csv_files"
    processes_no = 24
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
stop              = None
fasta             = "GAASMMGDVKESKMQITPETPGRIPVLNPFESPSDYSNLHEQTLASPSVFKSTKLPTPGKFRWSIDQLAVINPVEIDPEDIHRQALYLSHSRIDKDVEDKRQKAIEEFFTKDVIVPSPWTDHEGKQLSQCHSSKCTNINSDSPVGKKLTIHSEKSD"
folder2ptms       = get_folder2ptms(csvpath)
# for mz, intensity, q, path in islice(iter_scans(data_path), stop):
# mz, intensity, q, path = next(iter_scans(data_path))
# mz, intensity, q, path = next(filter(lambda x: len(folder2ptms[x[3].split("/")[-2]]) > 1, iter_scans(data_path)))
# it = iter_scans(data_path)
# mz, intensity, q, path = next(it)
# print(path)

def single_run(mz, intensity, q, path):
    exp, scan = path.split('/')[-2:]
    scan = int(scan)
    print("start: {exp}\tscan: {scan}".format(exp=exp, scan=scan))
    row = {"q":q, "exp":exp, "scan":scan}
    try:
        modifications = folder2ptms[exp]
        single_prec = len(modifications) in (0, 1)
        if single_prec:
            if modifications:
                modifications = modifications[0]
            else:
                modifications = {}
            M, timings = masstodon_single(mz, intensity, fasta, q, '',
                                          modifications     = modifications,
                                          isotopic_coverage = isotopic_coverage,
                                          min_prob          = min_prob,
                                          std_cnt           = std_cnt,
                                          orbitrap          = True,
                                          get_timings       = True,
                                          distance_charges  = 2)
        else:
            precursors = [{"name": "pBora-"+"_".join(map(str, mod.keys())),
                           "modifications": mod,
                           "q": q, 
                           "fasta": fasta,
                           "distance_charges": 2}
                          for mod in modifications]
            M, timings = masstodon_batch(mz, intensity, precursors,
                                          isotopic_coverage = isotopic_coverage,
                                          min_prob          = min_prob,
                                          std_cnt           = std_cnt,
                                          orbitrap          = True,
                                          get_timings       = True)
        local_dump_path = pjoin(dump_path, pjoin(exp), str(scan))
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
        if single_prec:
            for s in ('ETDorHTR', 'ETnoD_PTR_fragments', 'ETnoD_precursor', 'PTR_precursor'):
                row["cz_simple."+str(s)] = int(M.cz_simple.intensities[s])
                row["cz."+str(s)]        = int(M.cz.intensities[s])

        # estimates
        row['estimates'] = list(M.ome.iter_molecule_estimates())
        row['success'] = True
    except Exception as e:
        row['success'] = False
        raise e
    return row

ETD_data    = islice(filter(lambda x: "ETD" in x[3], iter_scans(data_path)), stop)
nonETD_data = islice(filter(lambda x: "ETD" not in x[3], iter_scans(data_path)), stop)


T0 = time()
with mp.Pool(processes_no) as p:
    stats = p.starmap(single_run, ETD_data)
T1 = time()
with open(pjoin(dump_path, "ETD_fit_stats.json"), "w") as f:
    json.dump(stats, f, indent=4)
print("Total fit time for all ETD spectra: {t}".format(t=T1-T0))


T2 = time()
with mp.Pool(processes_no) as p:
    stats = p.starmap(single_run, nonETD_data)
T3 = time()
with open(pjoin(dump_path, "nonETD_fit_stats.json"), "w") as f:
    json.dump(stats, f, indent=4)

print("Total fit time for all non-ETD spectra: {t}".format(t=T3-T2))
