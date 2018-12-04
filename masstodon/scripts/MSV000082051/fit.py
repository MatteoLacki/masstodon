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
experiment_specific_params = {}
Xs = []


X = "Oleic Acylation ETD 10 ms"
Xs.append(X)
experiment_specific_params[X] = dict(
    precursors = [dict(name = "ApoAI {} q={}".format(X,q),
                       fasta = ApoAI,
                       q = q,
                       distance_charges = distance_charges,
                       modifications = {133:{"C_alpha": ptms['oleic_acylation']}})
                  for q in charges])

X = "Carboxymethylation ETD 10 ms"
Xs.append(X)
experiment_specific_params[X] = dict(
    precursors = [dict(name = "ApoAI {} q={}".format(X,q),
                       fasta = ApoAI,
                       q = q, 
                       distance_charges = distance_charges,
                       modifications = {133:{"C_alpha": ptms['carboxymethylation']}})
                  for q in charges])

X = "Canonical ApoA-I ETD 5ms"
Xs.append(X)
experiment_specific_params[X] = dict(
    precursors = [dict(name = "ApoAI {} q={}".format(X,q),
                       fasta = ApoAI,
                       q = q, 
                       distance_charges = distance_charges)
                  for q in charges])

X = "ProtoApoA1 ETD 10ms"
Xs.append(X)
experiment_specific_params[X] = dict(
    precursors = [dict(name = "ProtoApoAI {} q={}".format(X,q),
                       fasta = ProtoApoAI,
                       q = q, 
                       distance_charges = distance_charges)
                  for q in charges])

X = "Arachidonic Acylation ETD 10ms"
Xs.append(X)
experiment_specific_params[X] = dict(
    precursors = [dict(name = "ApoAI {} q={}".format(X,q),
                       fasta = ApoAI,
                       q = q, 
                       distance_charges = distance_charges,
                       modifications = {88:{"C_alpha": ptms['arachidonic_acylation']}})
                  for q in charges])

X = "Arachidonic Acylation + Truncation ETD 7 ms"
Xs.append(X)
experiment_specific_params[X] = dict(
    precursors = [dict(name = "ApoAI {} q={}".format(X,q),
                       fasta = ApoAI[:242],
                       q = q, 
                       distance_charges = distance_charges,
                       modifications = {88:{"C_alpha": ptms['arachidonic_acylation']}})
                  for q in charges])

X = "Phosphorylation ETD 7 ms"
Xs.append(X)
experiment_specific_params[X] = dict(
    precursors = [dict(name = "ApoAI {} q={}".format(X,q),
                       fasta = ApoAI,
                       q = q, 
                       distance_charges = distance_charges,
                       modifications = {115:{"C_alpha": ptms['phosphorylation']}})
                  for q in charges])

X = "Oxidation ETD 5 ms"
Xs.append(X)
experiment_specific_params[X] = dict(
    precursors = [dict(name = "ApoAI {} q={}".format(X,q),
                       fasta = ApoAI,
                       q = q, 
                       distance_charges = distance_charges,
                       modifications = {86:{"C_alpha": ptms['oxidation']}})
                  for q in charges])

X = "Docohexaneoic Acylation ETD 10 ms"
Xs.append(X)
experiment_specific_params[X] = dict(
    precursors = [dict(name = "ApoAI {} q={}".format(X,q),
                       fasta = ApoAI,
                       q = q, 
                       distance_charges = distance_charges,
                       modifications = {88:{"C_alpha": ptms['docosahexanoic_acylation']}})
                  for q in charges])

X = "Palmitic Acylation + Truncation ETD 10 ms"
Xs.append(X)
experiment_specific_params[X] = dict(
    precursors = [dict(name = "ApoAI {} q={}".format(X,q),
                       fasta = ApoAI[:242],
                       q = q, 
                       distance_charges = distance_charges,
                       modifications = {88:{"C_alpha": ptms['palmitic_acylation']}})
                  for q in charges])

X = "ProtoApoA1 Oxidation ETciD 5ms 12ev"
Xs.append(X)
experiment_specific_params[X] = dict(
    precursors = [dict(name = "ProtoApoAI {} q={}".format(X,q),
                       fasta = ProtoApoAI,
                       q = q, 
                       distance_charges = distance_charges,
                       modifications = {92:{"C_alpha": ptms['oxidation']}})
                  for q in charges])

X = "Oleic Acylation + Truncation ETD 10 ms"
Xs.append(X)
experiment_specific_params[X] = dict(
    precursors = [dict(name = "ProtoApoAI {} q={}".format(X,q),
                       fasta = ProtoApoAI[:242],
                       q = q, 
                       distance_charges = distance_charges,
                       modifications = {88:{"C_alpha": ptms['oleic_acylation']}})
                  for q in charges])

X = "Dehydration ETD 10 ms"
Xs.append(X)
experiment_specific_params[X] = dict(
    precursors = [dict(name = "ApoAI {} q={}".format(X,q),
                       fasta = ApoAI,
                       q = q, 
                       distance_charges = distance_charges,
                       modifications = {PTM_position:{"C_alpha": ptms['dehydralation']}})
                  for q in charges for PTM_position in range(150, 167)])

X = "Truncation ETciD 5ms 12v"
Xs.append(X)
experiment_specific_params[X] = dict(
    precursors = [dict(name = "ApoAI {} q={}".format(X,q),
                       fasta = ApoAI[:242],
                       q = q, 
                       distance_charges = distance_charges)
                  for q in charges])

X = "Palmitic Acylation ETD 10 ms"
Xs.append(X)
experiment_specific_params[X] = dict(
    precursors = [dict(name = "ApoAI {} q={}".format(X,q),
                       fasta = ApoAI,
                       q = q, 
                       distance_charges = distance_charges,
                       modifications = {88:{"C_alpha": ptms['palmitic_acylation']}})
                  for q in charges])

X = "Glyacatio by Hexose ETD 7 ms"
Xs.append(X)
experiment_specific_params[X] = dict(
    precursors = [dict(name = "ApoAI {} q={}".format(X,q),
                       fasta = ApoAI,
                       q = q, 
                       distance_charges = distance_charges,
                       modifications = {131:{"C_alpha": ptms['glycation_by_hexose']}})
                  for q in charges])

for e in experiment_specific_params:
    experiment_specific_params[e]["mz"], experiment_specific_params[e]["intensity"] =\
        spectra[e]

