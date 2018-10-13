"""This needs revisiting. But the question is rather more about using
different error function altogether."""

from masstodon.masstodon import masstodon_single
from masstodon.data.substance_p import substance_p
from masstodon.data.constants import infinity


mz, intensity = substance_p['spectrum']

threshold         = 0.025
fasta             = 'RPKPQQFFGLM'
modifications     = {'11': {'C_carbo': {'H': 1, 'N': 1, 'O': -1}}}
# modifications = {}
q                 = 3
min_prob          = .8
isotopic_coverage = .999
std_cnt           = 3

m = masstodon_single(mz, intensity, fasta, q,
                     # min_mz_diff   = infinity,
                     modifications = modifications,
                     orbitrap      = False,
                     threshold     = threshold,
                     isotopic_coverage  = isotopic_coverage,
                     min_prob           = min_prob, 
                     std_cnt            = std_cnt,
                     include_zero_intensities = False)

fitted_dots = m.plotlygl("/Users/matteo/Projects/masstodon/masstodon/dump/deconvs_comparison/nnls")


n = masstodon_single(mz, intensity, fasta, q,
                     # min_mz_diff   = infinity,
                     modifications = modifications,
                     orbitrap      = False,
                     threshold     = threshold,
                     isotopic_coverage  = isotopic_coverage,
                     min_prob           = min_prob, 
                     std_cnt            = std_cnt,
                     include_zero_intensities = False,
                     deconvolution_method = 'quantile')

n.plotlygl("/Users/matteo/Projects/masstodon/masstodon/dump/deconvs_comparison/nnls")

sm = m.imperator.solutions[0]
sn = n.imperator.solutions[0]

sm.mean_mz
sn.mean_mz
sm.model.fitted()
sn.model.fitted()
sn.model.model['x']
sn.model.coef()
sn.fitted()[0]

import numpy as np

np.array((sn.model.X @ sn.model.coef())).flatten()

for sm, sn in zip(m.imperator.solutions, n.imperator.solutions):
    print(sm.model.fitted())
    print(sn.model.fitted())
    print()