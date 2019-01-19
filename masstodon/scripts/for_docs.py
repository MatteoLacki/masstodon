from pprint import pprint

from masstodon.data.constants import infinity
from masstodon.masstodon import masstodon_single
from masstodon.read.txt import spectrum_from_txt

path = "/Users/matteo/Projects/masstodon/data/Belgian/2013_05_01/SUBP/wave_height_1.5/FRL-010513-SUBP-WH 1,5-WV 4500.txt"

print("Opening spectrum.")
mz, i = spectrum_from_txt(path)

threshold         = "0.05 Th"
fasta             = 'RPKPQQFFGLM'
modifications     = {'11': {'C_carbo': {'H': 1, 'N': 1, 'O': -1}}}
q                 = 3
min_prob          = .8
isotopic_coverage = .999
std_cnt           = 3
timings           = True
include_zero_intensities = False


print("Running masstodon.")
m, t = masstodon_single(mz, i, fasta, q,
    min_mz_diff        = infinity,
    modifications      = modifications,
    orbitrap           = False,
    threshold          = threshold,
    isotopic_coverage  = isotopic_coverage,
    min_prob           = min_prob,
    std_cnt            = std_cnt,
    get_timings        = timings,
    include_zero_intensities = include_zero_intensities)

print("Save spectrum to where you started your python from.")
m.plotly("spectrum.html", shape='rectangles', show=False)

print("Getting stats on inintial and final number of nodes and edges in the deconvolution graph.")
pprint(m.ome.G_stats)

print("Errors: absolute deviations and relative distance between spectra.")
pprint(m.imperator.errors())

print("Estimates of intensities of reactions and fragmentations.")
pprint(m.cz_simple.intensities)
pprint(m.cz_simple.probabilities)

print("Timings.")
pprint(t)

print("That's pretty much it.")