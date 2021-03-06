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
from masstodon.masstodon import masstodon_single, masstodon_batch
from masstodon.data.ptms import ptms
from masstodon.precursor.ptm_filter import filter_ptm_assignments

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
threshold         = "5 ppm"
charges           = [35, 36, 37]
distance_charges  = 5
ApoAI             = "DEPPQSPWDRVKDLATVYVDVLKDSGRDYVSQFEGSALGKQLNLKLLDNWDSVTSTFSKLREQLGPVTQEFWDNLEKETEGLRQEMSKDLEEVKAKVQPYLDDFQKKWQEEMELYRQKVEPLRAELQEGARQKLHELQEKLSPLGEEMRDRARAHVDALRTHLAPYSDELRQRLAARLEALKENGGARLAEYHAKATEHLSTLSEKAKPALEDLRQGLLPVLESFKVSFLSALEEYTKKLNTQ"
ProtoApoAI        = "RHFWQQDEPPQSPWDRVKDLATVYVDVLKDSGRDYVSQFEGSALGKQLNLKLLDNWDSVTSTFSKLREQLGPVTQEFWDNLEKETEGLRQEMSKDLEEVKAKVQPYLDDFQKKWQEEMELYRQKVEPLRAELQEGARQKLHELQEKLSPLGEEMRDRARAHVDALRTHLAPYSDELRQRLAARLEALKENGGARLAEYHAKATEHLSTLSEKAKPALEDLRQGLLPVLESFKVSFLSALEEYTKKLNTQ"


experiments = {}
Xs = []
def insert_data(file,
                ptm_position=None,
                ptm=None,
                name="ApoAI {} q={}",
                fasta=ApoAI):
    if ptm_position is not None:
        assignments = filter_ptm_assignments(
            dict(name = name.format(file, q),
                 fasta = fasta,
                 q = q,
                 distance_charges = distance_charges,
                 modifications = {ptm_position:{"C_alpha":ptms[ptm]}})
            for q in charges)
        assignments = list(assignments)
    else:
        assignments = [
            dict(name = name.format(file, q),
                 fasta = fasta,
                 q = q,
                 distance_charges = distance_charges) 
            for q in charges]
    if len(assignments) > 0:
        Xs.append(file)
        experiments[file] = dict(precursors=assignments)
    else:
        print("There is no place for {} on {}.".format(ptm, file))

insert_data("Oleic Acylation ETD 10 ms",
            133, 'oleic_acylation')
insert_data("Carboxymethylation ETD 10 ms",
            133, 'carboxymethylation')
insert_data("Canonical ApoA-I ETD 5ms")
insert_data("ProtoApoA1 ETD 10ms",
    name="ProtoApoAI {} q={}",
    fasta=ProtoApoAI)
insert_data("Arachidonic Acylation ETD 10ms",
            88, 'arachidonic_acylation')
insert_data(
    "Arachidonic Acylation + Truncation ETD 7 ms",
    88, 'arachidonic_acylation',
    fasta = ApoAI[:242])
insert_data("Phosphorylation ETD 7 ms",
            115, 'phosphorylation')
insert_data("Oxidation ETD 5 ms",
            86, 'oxidation')
insert_data("Docohexaneoic Acylation ETD 10 ms",
            88, 'docosahexanoic_acylation')
insert_data(
    "Palmitic Acylation + Truncation ETD 10 ms",
    88, 'palmitic_acylation',
    fasta = ApoAI[:242])
insert_data("ProtoApoA1 Oxidation ETciD 5ms 12ev",
    92, 'oxidation',
    "ProtoApoAI {} q={}",
    ProtoApoAI)
insert_data("Oleic Acylation + Truncation ETD 10 ms",
    88, 'oleic_acylation',
    "ProtoApoAI {} q={}",
    ProtoApoAI[:242])
insert_data(
    "Truncation ETciD 5ms 12v",
    fasta = ApoAI[:242])
insert_data("Palmitic Acylation ETD 10 ms",
            88, 'palmitic_acylation')
insert_data("Glyacatio by Hexose ETD 7 ms",
            131, 'glycation_by_hexose')

# X = "Dehydration ETD 10 ms"
# Xs.append(X)
# experiments[X] = dict(
#     precursors=[dict(name = "ApoAI {} q={}".format(X,q),
#                      fasta = ApoAI,
#                      q = q, 
#                      distance_charges = distance_charges,
#                      modifications = {PTM_position:{"C_alpha": ptms['dehydralation']}})
#                   for q in charges for PTM_position in range(150, 167)])
# experiments[X]['precursors'] = list(filter_ptm_assignments(experiments[X]['precursors']))

for e in experiments:
    experiments[e]["mz"], experiments[e]["intensity"] =\
        spectra[e]

def iter_data(experiments):
    for exp, d in experiments.items():
        yield d['mz'], d['intensity'], exp, d['precursors']

# mz, intensity, exp, precursors = next(iter_data(experiments))
def single_run(mz, intensity, exp, precursors, verbose=True):
    if verbose:
        print("Running {}.".format(exp))
    out_path = pjoin(out_folder, exp.replace(" ", "_"))
    row = {"exp":exp}
    try:
        M, timings = masstodon_batch(
          mz, intensity, precursors,
          isotopic_coverage = isotopic_coverage,
          min_prob          = min_prob,
          std_cnt           = std_cnt,
          orbitrap          = orbitrap,
          threshold         = threshold,
          get_timings       = get_timings)
        mkdir(out_path, exist_ok=True)
        M.dump(out_path, indent=4)
        M.write(out_path)
        M.plotly(pjoin(out_path, "spectrum.html"),
                 show = False)
        with open(pjoin(out_path, 'timings.json'), 'w') as h:
            json.dump(timings, h, indent=4)
        row.update(M.imperator.errors())
        row.update({"t_"+str(n): T for n,T in timings})
        row['estimates'] = list(M.ome.iter_molecule_estimates())
        row['success'] = True
        if verbose:
            print('Finished with {}'.format(exp))
    except Exception as e:
        row['success'] = False
        raise e
        if verbose:
            print("Error in experiment {}.".format(exp))
    return row

T0 = time()
with mp.Pool(processes_no) as p:
    stats = p.starmap(single_run, iter_data(experiments))
T1 = time()
with open(pjoin(out_folder, "fit_stats.json"), "w") as f:
    json.dump(stats, f, indent=4)
print("Total fit time for all ETD spectra: {t}".format(t=T1-T0))

