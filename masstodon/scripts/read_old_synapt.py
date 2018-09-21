%load_ext autoreload
%autoreload 2
%load_ext line_profiler

from masstodon.read.txt import spectrum_from_txt
from masstodon.masstodon import masstodon_single

datapath = "/Users/matteo/Projects/masstodon/data/Belgian/2013_05_01/SUBP/wave_height_1.5/FRL-010513-SUBP-WH 1,5-WV 10.txt"
mz, intensity = spectrum_from_txt(datapath)

from masstodon.data.substance_p import substance_p

mz, intensity = substance_p['spectrum']

threshold = 0.05
fasta             = 'RPKPQQFFGLM'
modifications     = {'11': {'C_carbo': {'H': 1, 'N': 1, 'O': -1}}}
q                 = 3
min_prob          = .8
isotopic_coverage = .999
std_cnt           = 3

m = masstodon_single(mz, intensity, fasta, q, 
                     modifications = modifications,
                     orbitrap = False,
                     threshold=threshold,
                     isotopic_coverage = isotopic_coverage,
                     min_prob = min_prob, 
                     std_cnt  = std_cnt)

m.plotly("dump")


