%load_ext autoreload
%autoreload 2
%load_ext line_profiler

from collections import Counter
from os.path import join as pjoin
from sys import platform

from masstodon.masstodon import masstodon_single, Masstodon
from masstodon.read.npy  import spectrum_from_npy
from masstodon.scripts.PXD001845.get_fasta import read_fasta
from masstodon.scripts.PXD001845.hecks_custom_ptms import modify_fasta
from masstodon.scripts.PXD001845.csv2fasta import csv2fasta

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

scan = 10
folder = "20141202_AMB_Bora_10x_40MeOH_1FA_OT_120k_10uscans_920_EThcD_8ms_15CE_19precZ"

# folder = "20141202_AMB_pBora_PLK_10x_40MeOH_1FA_OT_120k_10uscans_928_ETciD_8ms_15SA_19precZ"
# scan = 1
# folder = "20141202_AMB_pBora_PLK_10x_40MeOH_1FA_OT_120k_10uscans_928_ETciD_8ms_15SA_19precZ"
# scan   = 1

full_data_path = pjoin(
    data_path, 
    "20141202_AMB_Bora_10x_40MeOH_1FA_OT_120k_10uscans_920_EThcD_8ms_15CE_19precZ/",
    str(scan))

mz, intensity = spectrum_from_npy(full_data_path)
# from masstodon.plot.spectrum import plot_spectrum
# plot_spectrum(mz, intensity)

min_prob          = .8
isotopic_coverage = .999
std_cnt           = 3
q                 = 19

matching_csv = [f for f in csv2fasta if folder in f]
fasta = csv2fasta[matching_csv[0]]
fasta, modifications = modify_fasta(fasta)


from masstodon.masstodon import Masstodon
from masstodon.spectrum.orbitrap import OrbitrapSpectrum
from masstodon.spectrum.cluster import Bitonic, Groups
from masstodon.stats.gaussian import mean, sd, skewness
from masstodon.spectrum.peak_clustering import bitonic_clustering, bitonic_iterator, fix_local_clustering

import warnings
import numpy as np
warnings.filterwarnings('error')

# m = Masstodon()
# m.set_spectrum(mz, intensity, orbitrap=True)

spec = OrbitrapSpectrum(mz, intensity)
# spec.plot()
# spec.bitonic_clustering()
# bc.fit(spec.mz, spec.intensity)
# bc.get_stats()

bc = Bitonic()
bc.x = spec.mz
bc.w = spec.intensity
spec.bc = bc

# this is how bc looks without correction:
# bc.clusters = np.array(list(bitonic_iterator(spec.mz, spec.intensity)))

# this is how it looks like corrected
bc.clusters = bitonic_clustering(bc.x, bc.w, .15, .2)
# spec.plot(clusters='bitonic')

from masstodon.models.polynomial import polynomial
from masstodon.stats.descriptive import mad_denoising

bc.fit_diff_model(model = polynomial,
                  degree= 3,
                  denoise_refit = {'denoiser': mad_denoising,
                                   'std_cnt' : 100})
bc.diff_model.plot()
bc.diff_model.coef()

import matplotlib.pyplot as plt
Z = np.array([(np.median(lmz), len(lmz)) for lmz, lintensity in bc])

plt.scatter(Z[:,0], Z[:,1])

plt.plot(np.arange(Z.shape[0]), Z[:,1])
plt.show()

Counter(Z[:,1])




# len(lefts)
# len(bc.clusters)


# diffs[0] == diffs[-1]

# lefts, diffs = bc.left_ends_and_diffs()
# np.percentile(diffs, 25)


# i = 0
# for local_mz, local_intensity in bc:
#     try:
#         with warnings.catch_warnings():
#             mean_mz = mean(local_mz, local_intensity)
#     except RuntimeWarning as e:
#         i += 1
#         print(local_mz)
#         print(local_intensity)
#         print()
#     # sd_mz         = sd(local_mz, local_intensity, mean_mz)
#     # skewnesses_mz = skewness(local_mz, local_intensity, mean_mz, sd_mz)
#     # min_mz        = min(local_mz)
#     # max_mz        = max(local_mz)
#     # for o, v in zip(O, (min_mz, max_mz, mean_mz, sd_mz, skewnesses_mz,\
#     #                     len(local_mz), sum(local_intensity))): 
#     #     o.append(v)
# # 61 cases of poor fitting.

# # testing bitonic clustering

# x = np.array([1,2,3,4,5,6,7])
# y = np.array([1,3,2,1,4,5,2])
# bitonic_clustering(x, y, 2)

# x = np.arange(1,8)
# y = np.array([1,3,2,6,4,5,2])
# bitonic_clustering(x, y, 2)

# x = np.array([1,2,3,5,6,7,8])
# y = np.array([1,3,2,1,4,5,2])
# bitonic_clustering(x, y, 3)
# list(bitonic_iterator(x, y, 3))






