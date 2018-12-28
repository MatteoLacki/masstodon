from itertools import islice
import json
import multiprocessing as mp
from os import listdir as ls, makedirs as mkdir
from os.path import join as pjoin
import re
from sys import platform
from time import time

from masstodon.data.ptms import ptms
from masstodon.data.constants              import infinity
from masstodon.deconvolve.divide_ed_impera import imperator, load_imperator
from masstodon.estimates_matcher.cz        import CzMatch
from masstodon.estimates_matcher.cz_simple import SimpleCzMatch
from masstodon.isotopes.calculator         import isotope_calculator
from masstodon.precursor.precursor         import precursor
from masstodon.preprocessing.filters       import filter_subspectra_molecules
from masstodon.read.npy                    import spectrum_from_npy
from masstodon.spectrum.spectrum           import spectrum
from masstodon.ome.ome                     import ome

class WrongSystem(Exception):
    pass

if platform == "darwin": # latte on skimmed soya milk anyone?
    path = "/Users/matteo/Projects/masstodon/data/MSV000082051/mzmls/"
    out_folder = "/Users/matteo/Projects/masstodon/data/MSV000082051/res"
    upper_limit_of_masstodon_runs = 4
    processes_no = 4
elif platform == "linux": # death metal and long dirty hair, fuck yeah!
    stem_path = "/mnt/disk/masstodon/data/MSV000082051"
    path = pjoin(stem_path, "mzmls")
    csvpath = pjoin(stem_path, "csvs")
    out_folder = pjoin(stem_path, "res")
    upper_limit_of_masstodon_runs = None
    processes_no = 24
else:
    raise WrongSystem("Path not specified correctly.")

files = [x for x in ls(path) if ".mzML" in x]
spectra = {f.replace(".mzML",""): (s['m/z array'], s['intensity array'])
           for f in files for s in read_mzml(pjoin(path,f)) if s['ms level'] == 2}
files = list(spectra)

min_prob          = .8
isotopic_coverage = .999
std_cnt           = 3
get_timings       = True
orbitrap          = True
charges           = [35, 36, 37]
distance_charges  = 5
ApoAI             = "DEPPQSPWDRVKDLATVYVDVLKDSGRDYVSQFEGSALGKQLNLKLLDNWDSVTSTFSKLREQLGPVTQEFWDNLEKETEGLRQEMSKDLEEVKAKVQPYLDDFQKKWQEEMELYRQKVEPLRAELQEGARQKLHELQEKLSPLGEEMRDRARAHVDALRTHLAPYSDELRQRLAARLEALKENGGARLAEYHAKATEHLSTLSEKAKPALEDLRQGLLPVLESFKVSFLSALEEYTKKLNTQ"
ProtoApoAI        = "RHFWQQDEPPQSPWDRVKDLATVYVDVLKDSGRDYVSQFEGSALGKQLNLKLLDNWDSVTSTFSKLREQLGPVTQEFWDNLEKETEGLRQEMSKDLEEVKAKVQPYLDDFQKKWQEEMELYRQKVEPLRAELQEGARQKLHELQEKLSPLGEEMRDRARAHVDALRTHLAPYSDELRQRLAARLEALKENGGARLAEYHAKATEHLSTLSEKAKPALEDLRQGLLPVLESFKVSFLSALEEYTKKLNTQ"
experiments = {}
Xs = []

# The biggest copy-paste ever!!!
X = "Oleic Acylation ETD 10 ms"
Xs.append(X)
experiments[X] = dict(
    precursors = [dict(name = "ApoAI {} q={}".format(X,q),
                       fasta = ApoAI,
                       q = q,
                       distance_charges = distance_charges,
                       modifications = {133:{"C_alpha": ptms['oleic_acylation']}})
                  for q in charges])
mz, intensity = spectra[X]
S = spectrum(mz, intensity, 1.1, orbitrap)
# S.plot()
subspectra = S.get_min_mz_diff_subspectra()

mz_digits = S.get_smallest_diff_digits()
iso_calc = isotope_calculator(digits = mz_digits)
std_cnt = 3



P = precursor(iso_calc=iso_calc, **experiments[X]['precursors'][2])
mols = list(P.molecules())
len(mols)
O = ome(iso_calc, experiments[X]['precursors'])
good_mols, good_subspectra = O.filter_by_deviations(subspectra)
O.G_stats

