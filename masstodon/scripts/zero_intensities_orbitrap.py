%load_ext autoreload
%autoreload 2
%load_ext line_profiler

from time import time
from collections import Counter

from masstodon.masstodon import Masstodon, masstodon_batch, masstodon_single, load_masstodon
from masstodon.ome.ome   import Ome
from masstodon.read.npy  import spectrum_from_npy

data_path = '/Users/matteo/Projects/masstodon/data/PXD001845/numpy_files/20141202_AMB_pBora_PLK_10x_40MeOH_1FA_OT_120k_10uscans_928_ETciD_8ms_15SA_19precZ/1'
# data_path = "/Users/matteo/Projects/masstodon/data/PXD001845/numpy_files/20141202_AMB_Bora_10x_40MeOH_1FA_OT_120k_10uscans_920_ETD_4ms_19precZ/6"

mz, intensity = spectrum_from_npy(data_path)
fasta             = "GAASMMGDVKESKMQITPETPGRIPVLNPFESPSDYSNLHEQTLASPSVFKSTKLPTPGKFRWSIDQLAVINPVEIDPEDIHRQALYLSHSRIDKDVEDKRQKAIEEFFTKDVIVPSPWTDHEGKQLSQCHSSKCTNINSDSPVGKKLTIHSEKSD"
charge            = 24
min_prob          = .8
isotopic_coverage = .999
std_cnt           = 3

m, tm = masstodon_single(mz, intensity, fasta, charge, '',
                         orbitrap = True,
                         isotopic_coverage = isotopic_coverage,
                         min_prob = min_prob, 
                         std_cnt  = std_cnt,
                         get_timings = True,
                         include_zero_intensities = True)
m.plotlygl("dump/mit_zeros")

n, tn = masstodon_single(mz, intensity, fasta, charge, '',
                     orbitrap = True,
                     isotopic_coverage = isotopic_coverage,
                     min_prob = min_prob, 
                     std_cnt  = std_cnt,
                     get_timings = True,
                     include_zero_intensities = False)
n.plotlygl("dump/ohne_zeros")

