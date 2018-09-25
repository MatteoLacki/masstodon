from masstodon.masstodon import masstodon_single
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
                     min_mz_diff   = infinity,
                     modifications = modifications,
                     orbitrap      = False,
                     threshold     = threshold,
                     isotopic_coverage = isotopic_coverage,
                     min_prob = min_prob, 
                     std_cnt  = std_cnt)

m.dump("")
m.plotlygl("", show=False)