%load_ext autoreload
%autoreload 2
from itertools import islice
import json
import multiprocessing as mp
from os import listdir as ls, makedirs as mkdir
from os.path import join as pjoin
import re
from sys import platform
from time import time

from masstodon.read.mzml import read_mzml
from masstodon.plot.spectrum import plot_spectrum
from masstodon.scripts.PXD001845.parse_csvs_with_precursors2 import get_folder2ptms
from masstodon.masstodon import masstodon_single, masstodon_batch

class WrongSystem(Exception):
    pass

if platform == "darwin": # latte on skimmed soya milk anyone?
    path = "/Users/matteo/Projects/masstodon/data/PXD001845/mzmls/"
    csvpath = "/Users/matteo/Projects/masstodon/data/PXD001845/csvs"
    out_folder = "/Users/matteo/Projects/masstodon/data/PXD001845/res"
    upper_limit_of_masstodon_runs = 4
    processes_no = 4
elif platform == "linux": # death metal and long dirty hair, fuck yeah!
    stem_path = "/mnt/disk/masstodon/data/PXD001845"
    path = pjoin(stem_path, "mzmls")
    csvpath = pjoin(stem_path, "csvs")
    upper_limit_of_masstodon_runs = None
    processes_no = 24
else:
    raise WrongSystem("Path not specified correctly.")

files = [x for x in ls(path) if ".mzML" in x]
ETD = [f for f in files if "ETD" in f] # independent, multicore execution of masstodons
min_prob          = .8
isotopic_coverage = .999
std_cnt           = 3
stop              = None
fasta             = "GAASMMGDVKESKMQITPETPGRIPVLNPFESPSDYSNLHEQTLASPSVFKSTKLPTPGKFRWSIDQLAVINPVEIDPEDIHRQALYLSHSRIDKDVEDKRQKAIEEFFTKDVIVPSPWTDHEGKQLSQCHSSKCTNINSDSPVGKKLTIHSEKSD"
folder2ptms       = get_folder2ptms(csvpath)
get_timings       = True
orbitrap          = True
# this has to be enhanced for data other than ETD:
get_q = lambda f: int(re.sub('[^0-9]','', f.split("_")[-1]))

def iter_data(out_folder, files=ETD):
    for f in files:
        q = get_q(f)
        exp = f.replace('.mzML','')
        mods = folder2ptms.get(exp,[])
        for s in read_mzml(pjoin(path,f)):
            scan_no   = int(s['id'].split("=")[1])
            mz        = s['m/z array']
            intensity = s['intensity array']
            yield mz, intensity, exp, scan_no, q, fasta, mods

# mz, intensity, exp, scan_no, q, fasta, mods = next(iter_data(out_folder))
def single_run(mz, intensity, exp, scan_no, q, fasta, mods, verbose=True):
    if verbose:
        print(exp, scan_no)
    out_path = pjoin(out_folder, exp, str(scan_no))
    row = {"q":q, "exp":exp, "scan":scan_no}
    try:
        modifications = folder2ptms.get(exp, [])
        single_prec = len(modifications) in (0, 1)
        if single_prec: # one PTM == no PTM-free
            modifications = modifications[0] if modifications else {}
            M, timings = masstodon_single(mz, intensity, fasta, q, '',
                                          modifications     = modifications,
                                          isotopic_coverage = isotopic_coverage,
                                          min_prob          = min_prob,
                                          std_cnt           = std_cnt,
                                          orbitrap          = orbitrap,
                                          get_timings       = get_timings,
                                          distance_charges  = 2)
        else:
            precursors = [{"name": "pBora-"+"_".join(map(str, mod.keys())),
                           "modifications": mod,
                           "q": q, 
                           "fasta": fasta,
                           "distance_charges": 2}
                          for mod in modifications] # phosphorylation == pBora
            M, timings = masstodon_batch(mz, intensity, precursors,
                                          isotopic_coverage = isotopic_coverage,
                                          min_prob          = min_prob,
                                          std_cnt           = std_cnt,
                                          orbitrap          = orbitrap,
                                          get_timings       = get_timings)
        mkdir(out_path, exist_ok=True)
        M.dump(out_path, indent=4)
        M.write(out_path)
        M.plotly(pjoin(out_path, "spectrum.html"), show=False)
        with open(pjoin(out_path, 'timings.json'), 'w') as h:
            json.dump(timings, h, indent=4)
        row.update(M.imperator.errors())
        row.update({"t_"+str(n): T for n,T in timings})
        if single_prec: # intensities
            for s in ('ETDorHTR', 'ETnoD_PTR_fragments', 'ETnoD_precursor', 'PTR_precursor'):
                row["cz_simple."+str(s)] = int(M.cz_simple.intensities[s])
                row["cz."+str(s)] = int(M.cz.intensities[s])
        row['estimates'] = list(M.ome.iter_molecule_estimates())
        row['success'] = True
    except Exception as e:
        row['success'] = False
        raise e
    return row

data = islice(iter_data(out_folder, files=ETD), upper_limit_of_masstodon_runs)
T0 = time()
with mp.Pool(processes_no) as p:
    stats = p.starmap(single_run, data)
T1 = time()
with open(pjoin(out_folder, "fit_stats.json"), "w") as f:
    json.dump(stats, f, indent=4)
print("Total fit time for all ETD spectra: {t}".format(t=T1-T0))

