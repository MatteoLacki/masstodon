# %load_ext autoreload
# %autoreload 2
# %load_ext line_profiler

from   collections  import  defaultdict, namedtuple, Counter
import numpy        as      np
import networkx     as      nx
try:
    import matplotlib.pyplot as plt
except RuntimeError:
    pass
from   time                 import  time
from   math                 import  log10, floor

from masstodon.readers.from_npy            import spectrum_from_npy
from masstodon.precursor.precursor         import precursor
from masstodon.isotopes                    import isotope_calculator
from masstodon.spectrum.spectrum           import spectrum
from masstodon.models.polynomial           import polynomial
from masstodon.preprocessing.filters       import filter_subspectra_molecules
from masstodon.deconvolve.divide_ed_impera import divide_ed_impera
from masstodon.estimates_matcher.cz_simple import SimpleCzMatch

if __name__ == '__main__':
    # generating subspectra
    data_path     = '/Users/matteo/Projects/review_masstodon/data/PXD001845/numpy_files/20141202_AMB_pBora_PLK_10x_40MeOH_1FA_OT_120k_10uscans_928_ETciD_8ms_15SA_19precZ/1'
    mz, intensity = spectrum_from_npy(data_path)

    spec = spectrum(mz, intensity)
    spec.bitonic_clustering()
    spec.min_mz_diff_clustering()

    # spec.plot()

    # spec.plot_mz_diffs()
    # spec.plot(clusters='bitonic')
    # spec.plot(clusters='min_mz_diff')
    # spec.bc.plot_sd()

    subspectra = list(spec.iter_min_mz_diff_subspectra())
    mz_digits  = spec.bc.get_smallest_diff_digits()
    iso_calc   = isotope_calculator(digits = mz_digits)

    # generating formulas
    fasta  = "GAASMMGDVKESKMQITPETPGRIPVLNPFESPSDYSNLHEQTLASPSVFKSTKLPTPGKFRWSIDQLAVINPVEIDPEDIHRQALYLSHSRIDKDVEDKRQKAIEEFFTKDVIVPSPWTDHEGKQLSQCHSSKCTNINSDSPVGKKLTIHSEKSD"
    charge = 24
    prec   = precursor(fasta, charge, name="", iso_calc=iso_calc)
    mols   = np.array(list(prec.molecules()))
    min_prob = .8
    isotopic_coverage = .999
    good_mols, good_subspectra = filter_subspectra_molecules(subspectra,
                                                             mols,
                                                             std_cnt = 3)

    t0 = time()
    imperator = divide_ed_impera(good_mols, spec.bc, min_prob, isotopic_coverage)
    fit_time = time() - t0

    # getting the biggest optimization problem available.
    a, i = max([ (s.XY()[0].shape[1], i) for i, s in enumerate(imperator.solutions)])
    X, Y = imperator.solutions[i].XY()
    X.shape

    from masstodon.stats.weighted_median import weighted_median


    np.all(X >= 0)
    X.shape, Y.shape
    from itertools import combinations



    def l1(a, b):
        return np.abs(a - b).sum()

    l1_distances = {}
    for i, j in combinations(range(X.shape[1]), 2):
        l1_distances[(i,j)] = l1(X[:,i], X[:,j])

    np.percentile(list(l1_distances.values()),
                  q = range(0,100, 10))

    min((v, k) for k, v in l1_distances.items())

    X[:,1]
    X[:,11]


    weighted_median()




    # imperator.plot()
    # imperator.plot_ccs()
    # imperator.plot_solutions(plt_style = 'ggplot')
    # imperator.solutions[10].plot(plt_style = 'fast')
    # imperator.solutions[10].plot_fancy(plt_style = 'ggplot')

