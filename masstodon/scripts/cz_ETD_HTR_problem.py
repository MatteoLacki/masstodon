"""What needs to be done here:

    We have to modify the pairing graph to discern in certain cases
    between ETD and HTR (it is possible).

    Right now, it does assume that HTR and ETD come together.
    Need to establish good criterion for matching ions and
    revisit the pairing procedure.
"""

%load_ext autoreload
%autoreload 2
%load_ext line_profiler

import pandas as pd
from os         import listdir, makedirs
from os.path    import join as pjoin, exists as pexists

from masstodon.read.txt     import spectrum_from_txt
from masstodon.masstodon    import masstodon_single
from masstodon.data.constants import infinity

common_path = "/Users/matteo/Projects/masstodon/data/Belgian/2013_05_01/SUBP"
folders = ("wave_velocity_300", "wave_height_1.5")

def iter_data(path, folders):
    """Iterate over existing spectra in different folders with common root at the path."""
    for f in folders:
        pf = pjoin(path, f)
        for cff in listdir(pf):
            cffp = pjoin(pf, cff)
            try:
                yield cff, spectrum_from_txt(cffp)
            except Exception as e:
                pass


dump_folder = "/Users/matteo/Projects/masstodon/dumps/belgian/synapt"

threshold         = 0.025
fasta             = 'RPKPQQFFGLM'
modifications     = {'11': {'C_carbo': {'H': 1, 'N': 1, 'O': -1}}}
q                 = 3
min_prob          = .8
isotopic_coverage = .999
std_cnt           = 3
timings           = True
include_zero_intensities = True

bad = []

def res_with_figs():
    for i, (cff, (mz, intensity)) in enumerate(iter_data(common_path, folders)):
        try:
            m, t = masstodon_single(mz, intensity, fasta, q,
                                 min_mz_diff        = infinity,
                                 modifications      = modifications,
                                 orbitrap           = False,
                                 threshold          = threshold,
                                 isotopic_coverage  = isotopic_coverage,
                                 min_prob           = min_prob, 
                                 std_cnt            = std_cnt,
                                 get_timings        = timings,
                                 include_zero_intensities = include_zero_intensities)
            df = pjoin(dump_folder, "{i}_".format(i=i) + cff)
            if not pexists(df):
                makedirs(df)
            m.dump(df)
            m.write(df)
            m.plotlygl(df, shape='rectangles', show=False)
            print(f'Finished with {cff}')
        except AssertionError:
            bad.append((cff, mz, intensity))

res_with_figs()
bad

# solving the issue of poor cz results.
from masstodon.masstodon import Masstodon
from masstodon.estimates_matcher.cz_simple import SimpleCzMatch
from masstodon.estimates_matcher.cz import CzMatch


cff, mz, intensity = bad[0]
m, t = masstodon_single(mz, intensity, fasta, q,
                     min_mz_diff        = infinity,
                     modifications      = modifications,
                     orbitrap           = False,
                     threshold          = threshold,
                     isotopic_coverage  = isotopic_coverage,
                     min_prob           = min_prob, 
                     std_cnt            = std_cnt,
                     get_timings        = timings,
                     include_zero_intensities = include_zero_intensities)

m = Masstodon()
m.set_spectrum(mz, intensity, infinity, False, threshold)
m.set_isotopic_calculator()
m.set_ome([{'fasta':fasta, 'q':q, 'modifications':modifications}])
m.divide_et_impera(min_prob, isotopic_coverage, include_zero_intensities)

prec = next(m.ome.sources())
list(prec.uncharged_molecules())
list(prec.molecules())

list(prec._protonate('p'))
list(prec._protonate('c'))
list(prec._protonate('z'))
from masstodon.precursor.precursor import precursor
P = precursor(fasta, 4, distance_charges=0)
list(P._protonate("p"))
list(P._protonate("c"))
list(P._protonate("z"))

sources = m.ome.sources()
prec = next(sources)
try:
    x = next(sources)
    raise AttributeError("You supplied too many precursors for the c/z analysis.")
except StopIteration:
    pass
for mol in m.ome.observables():
    mol.name = m.ome.G[mol][prec]['name']
    print(mol.name)
    mol.prec_fasta_len = len(prec.fasta)
prec = next(m.ome.sources())


m.ome.filter_by_estimated_intensity()
m.restrict_good_mols()

cz_simple = SimpleCzMatch(m.good_mols, prec.q)
cz = CzMatch(m.good_mols, prec.q)
cz.graph.nodes(data=True)

