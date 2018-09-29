%load_ext autoreload
%autoreload 2
%load_ext line_profiler

from os.path import join as pjoin
from sys import platform
import matplotlib.pyplot as plt

from masstodon.plot.spectrum import plot_spectrum
from masstodon.masstodon import masstodon_single, Masstodon
from masstodon.read.npy  import spectrum_from_npy
from masstodon.scripts.PXD001845.get_fasta import read_fasta
from masstodon.scripts.PXD001845.hecks_custom_ptms import modify_fasta
from masstodon.scripts.PXD001845.csv2fasta import csv2fasta
from masstodon.spectrum.orbitrap import OrbitrapSpectrum
from masstodon.models.polynomial import polynomial
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter


# data_path = "/Users/matteo/Projects/masstodon/data/PXD001845/numpy_files/20141202_AMB_pBora_AurA_10x_40MeOH_1FA_OT_60k_10uscans_924_EThcD_6ms_8CE_19precZ/15/"
# data_path = "/Users/matteo/Projects/masstodon/data/PXD001845/numpy_files/20141202_AMB_Bora_10x_40MeOH_1FA_OT_120k_10uscans_795_ETciD_4ms_15SA_22precZ/1/"
# data_path = "/Users/matteo/Projects/masstodon/data/PXD001845/numpy_files/20141202_AMB_Bora_10x_40MeOH_1FA_OT_120k_10uscans_728_ETD_4ms_24precZ/1/"
# data_path = "/Users/matteo/Projects/masstodon/data/PXD001845/numpy_files/20141202_AMB_Bora_10x_40MeOH_1FA_OT_120k_10uscans_728_ETD_4ms_24precZ/1/"

# czczmiel data
# data_path = "/Users/matteo/duch/20141202_AMB_Bora_10x_40MeOH_1FA_OT_120k_10uscans_728_ETD_4ms_24precZ/1/"
# data_path = "/home/matteo/masstodon/review_answer/numpy_files/20141202_AMB_Bora_10x_40MeOH_1FA_OT_120k_10uscans_728_ETD_4ms_24precZ/1/"

if platform == "darwin":
    # check if you have your latte on skimmed soya milk with you
    data_path = "/Users/matteo/Projects/masstodon/data/PXD001845/numpy_files/"
    dump_path = "/Users/matteo/Projects/masstodon/dumps/tests/"
elif platform == "linux":
    # check if you have long dirty hair
    data_path = "/home/matteo/masstodon/review_answer/numpy_files/"
    dump_path = "/mnt/disk/masstodon/dumps/tests/"

# "20141202_AMB_Bora_10x_40MeOH_1FA_OT_120k_10uscans_920_EThcD_8ms_15CE_19precZ/4/"
# folder = "20141202_AMB_Bora_10x_40MeOH_1FA_OT_120k_10uscans_920_EThcD_8ms_15CE_19precZ"
# scan = 1
# folder = "20141202_AMB_Bora_10x_40MeOH_1FA_OT_120k_10uscans_920_EThcD_8ms_15CE_19precZ"
# scan = 2

folder = "20141202_AMB_Bora_10x_40MeOH_1FA_OT_120k_10uscans_920_EThcD_8ms_15CE_19precZ"
scan = 2

full_data_path = pjoin(data_path, folder, str(scan))
mz, intensity = spectrum_from_npy(full_data_path)


min_prob          = .8
isotopic_coverage = .999
std_cnt           = 3
q                 = 19
matching_csv = [f for f in csv2fasta if folder in f]
fasta = csv2fasta[matching_csv[0]]
fasta, modifications = modify_fasta(fasta)
m,t = masstodon_single(mz, intensity, fasta, q, '',
                     modifications = modifications,
                     orbitrap = True,
                     isotopic_coverage = isotopic_coverage,
                     min_prob = min_prob, 
                     std_cnt  = std_cnt,
                     include_zero_intensities = False,
                     get_timings = True)
m.dump(dump_path)
m.plotlygl(dump_path)


# from masstodon.plot.spectrum import plot_spectrum
# plot_spectrum(mz, intensity)

spec = OrbitrapSpectrum(mz, intensity, drop_zeros=False)
# spec.plot(show=False)
# plt.scatter(mz[intensity==0], intensity[intensity==0])
# plt.show()

spec.bitonic_clustering()

# spec.plot(clusters='bitonic')

spec.bc.fit_diff_model(polynomial, degree = 2)
spec.bc.diff_model.plot()

x, y = spec.bc.diff_model.x, spec.bc.diff_model.y
res = spec.bc.diff_model.res()
order_lower  = np.median(res)/10 >= res
order_higher = np.median(res)*10 <= res

from masstodon.spectrum.lightweight import lightweight_spectrum

spec.bc.get_stats()
spec.bc.clusters

len(spec.bc)

np.arange(len(spec.bc.groups.max_mz))
lightweight_spectrum(spec.bc.groups.min_mz, spec.bc.groups.max_mz)

sum(order_lower)
# place this in a subclass of bitonic clustering

def iter_clusters_under(x):
    X = x[order_lower]
    i = 0
    Xi = X[i]
    for lmz, lintensity in spec.bc:
        if Xi in lmz:
            yield lmz, lintensity, i, Xi
            i += 1
            Xi = X[i]

under_clusters = list(iter_clusters_under(x))
len(under_clusters)




