"""This needs revisiting. But the question is rather more about using
different error function altogether."""

%load_ext autoreload
%autoreload 2

from masstodon.masstodon import masstodon_single, masstodon_batch
from masstodon.data.substance_p import substance_p
from masstodon.data.constants import infinity

mz, intensity = substance_p['spectrum']
threshold         = 0.025
fasta             = 'RPKPQQFFGLM'
modifications     = {'11': {'C_carbo': {'H': 1, 'N': 1, 'O': -1}}}
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
# prec = next(m.ome.sources())
# for n in m.ome.observables():
#     print(m.ome.G[prec][n]["name"])
#     print(n.formula.tex_with_charges(ce=True))
#     print()

m.dump("")
m.plotly(show=True, path='/Users/matteo/Projects/masstodon/substance_p.html')

precursors = [{"name": "substanceP",
               "modifications": modifications,
               "q": q, 
               "fasta": fasta,
               "distance_charges": 2}]
mb = masstodon_batch(
    mz,
    intensity,
    precursors = precursors)
