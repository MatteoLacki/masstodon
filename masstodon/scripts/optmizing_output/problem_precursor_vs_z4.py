%load_ext autoreload
%autoreload 2

from masstodon.data.substanceP_wh15_wv400 import mz, intensity, fasta, modifications, z
from masstodon.plot.spectrum    import plot_spectrum
from masstodon.masstodon        import masstodon_single
from masstodon.data.constants   import infinity

m = masstodon_single(mz, intensity, fasta, z,
                     modifications      = modifications,
                     orbitrap           = False,
                     threshold          = "0.02Da",
                     isotopic_coverage  = .999,
                     min_prob           = .8, 
                     std_cnt            = 3,
                     include_zero_intensities = False,
                     deconvolution_method = "nnls")


list(m.ome.iter_estimates_no_g(header=True))
list(m.ome.iter_estimates(header=True))

