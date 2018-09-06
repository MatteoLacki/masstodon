%load_ext autoreload
%autoreload 2
%load_ext line_profiler

from   collections          import  defaultdict, namedtuple, Counter
import numpy                as      np
import networkx             as      nx
import matplotlib.pyplot    as      plt
from   time                 import  time
from   math                 import  log10, floor

from masstodon.readers.from_npy               import spectrum_from_npy
from masstodon.precursor.precursor            import precursor
from masstodon.isotopes                       import isotope_calculator
from masstodon.spectrum.spectrum              import spectrum
from masstodon.models.polynomial              import polynomial
from masstodon.preprocessing.filters          import filter_subspectra_molecules
from masstodon.deconvolve.divide_ed_impera    import divide_ed_impera
from masstodon.estimates_matcher.cz_match_simple import SimpleCzMatch

# generating subspectra
data_path     = '/Users/matteo/Projects/review_masstodon/data/PXD001845/numpy_files/20141202_AMB_pBora_PLK_10x_40MeOH_1FA_OT_120k_10uscans_928_ETciD_8ms_15SA_19precZ/1'
mz, intensity = spectrum_from_npy(data_path)

spec = spectrum(mz, intensity)
spec.bitonic_clustering()
spec.min_mz_diff_clustering()

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


# imperator.plot()
# imperator.plot_ccs()
# imperator.plot_solutions()
# imperator.solutions[10].plot()
# imperator.solutions[10].plot_fancy()

