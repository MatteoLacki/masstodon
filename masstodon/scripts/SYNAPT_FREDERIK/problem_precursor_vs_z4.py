%load_ext autoreload
%autoreload 2

import numpy    as np
import pandas   as pd
from os         import listdir, makedirs
from os.path    import join as pjoin, exists as pexists

from masstodon.data.substanceP_wh15_wv400 import mz, intensity, fasta, modifications, z
from masstodon.plot.spectrum    import plot_spectrum
from masstodon.masstodon        import masstodon_single
from masstodon.data.constants   import infinity

m = masstodon_single(mz, intensity, fasta, z,
                     modifications      = modifications,
                     orbitrap           = False,
                     threshold          = 0.02,
                     isotopic_coverage  = .999,
                     min_prob           = .8, 
                     std_cnt            = 3,
                     include_zero_intensities = False,
                     deconvolution_method = "nnls")

m.plotly("/Users/matteo/Projects/masstodon/tests/spectrum.html", shape='rectangles', show=True)
m.iso_calc._isotope_DB
m.iso_calc(formula='C63H100N18O13S1', q=3, g=-2).mz

# what are the bins like???
bin_sizes = m.imperator.groups.max_mz - m.imperator.groups.min_mz
np.quantile(bin_sizes, q = np.array([0,.5,1]))
bin_sizes[bin_sizes < .049]


n = masstodon_single(mz, intensity, fasta, z,
                     modifications      = modifications,
                     orbitrap           = False,
                     threshold          = .05,
                     isotopic_coverage  = .999,
                     min_prob           = .8, 
                     std_cnt            = 3,
                     include_zero_intensities = False,
                     deconvolution_method = "nnls")

m2 = masstodon_single(mz, intensity, fasta, z,
                     modifications      = modifications,
                     orbitrap           = False,
                     threshold          = 0.025,
                     isotopic_coverage  = .999,
                     min_prob           = .8, 
                     std_cnt            = 3,
                     include_zero_intensities = False,
                     deconvolution_method = "quantile")

def get_z4(masstodon):
    for M in masstodon.ome.observables():
        if M.name == 'z4':
            z4_intensity = M.intensity
            break
    return z4_intensity

def write_estimates(masstodon):
    for M in masstodon.ome.observables():
        print(M.name, M.intensity)

write_estimates(m)
m_z4 = get_z4(m)
get_z4(m2)
n_z4 = get_z4(n)

m.imperator.solutions[1].cc.nodes
n.imperator.solutions[1].cc.nodes
# repeating fitting like m from the very beginning:

from masstodon.spectrum.spectrum           import spectrum
from masstodon.isotopes.calculator         import isotope_calculator
from masstodon.ome.ome                     import ome
from masstodon.deconvolve.divide_ed_impera import imperator

precursors = [{"fasta":               fasta,
               "q":                   z,
               "name":                "subP",
               "modifications":       modifications,
               "blocked_fragments":   ["c0"],
               "block_prolines":      True,
               "distance_charges":    5}]

spec_m = spectrum(mz, intensity, 1.1, False, .02, sort=True, drop_duplicates=True, drop_zeros=True)
subspectra_m = spec_m.get_min_mz_diff_subspectra()


mz_digits_m = spec_m.get_smallest_diff_digits()
ls_m        = spec_m.get_lightweight_spectrum()
groups_m    = spec_m.get_groups()
iso_calc_m  = isotope_calculator(digits=mz_digits_m)
ome_m       = ome(iso_calc_m, precursors)
good_mols_m, subspectra_within_mz_deviations_m = ome_m.filter_by_deviations(subspectra_m, 3)

for M in ome_m.G:
    print(M, ome_m.G.node[M])

imperator_m = imperator(good_mols_m,
                        groups_m,
                        ls_m,
                        .8,
                        .999,
                        "nnls",
                        False)

from    networkx.linalg.attrmatrix  import attr_matrix

main_clust = imperator_m.solutions[0]
main_clust.cc.nodes
main_clust.idx
main_clust.X
main_clust.model.coef()
list(main_clust.iter_estimates())

attr_matrix(main_clust.cc, edge_attr='prob')

total_intensities = groups_m.intensity
min_mz = groups_m.min_mz
max_mz = groups_m.max_mz
mean_mz = groups_m.mean_mz

include_000 = False
cc          = main_clust.cc
cc.edges(data=True)
cc.nodes

X, ordering = attr_matrix(cc, edge_attr='prob')
ordering    = np.array(ordering)
mol_columns = ordering < 0
peak_rows   = ordering >= 0
X           = X[:,mol_columns][peak_rows,:]
idx         = ordering[peak_rows]
            ordering[mol_columns]
total_intensities = total_intensities[idx]
mz_s        = min_mz[idx]
mz_e        = max_mz[idx]
mean_mz     = mean_mz[idx]
Y           = total_intensities
if include_000:
    Y1 = np.concatenate((Y, np.zeros(X.shape[1])))
    x  = 1.0 - np.array(X.sum(axis=0)).flatten()
    X1 = np.concatenate((X, np.diag(x)))
    model = fit_model(X1, Y1)
else:
    model = fit_model(X, Y)


mol_idx = ordering[mol_columns]
coefs   = main_clust.model.coef()

list(zip(mol_idx, coefs))






# imperator_m.mols[2]
# imperator_m.mols[2]
# imperator_m.mols[19]
# imperator_m.mols[40]
# list(ome_m.iter_nodes(False))
# spec_n = spectrum(mz, intensity, 1.1, False, .05, sort=True, drop_duplicates=True, drop_zeros=True)
# mz_digits_n = spec_n.get_smallest_diff_digits()
# ls_n        = spec_n.get_lightweight_spectrum()
# groups_n    = spec_n.get_groups()
# iso_calc_n  = isotope_calculator(digits=mz_digits_n)

