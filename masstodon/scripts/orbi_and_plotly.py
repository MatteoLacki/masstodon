%load_ext autoreload
%autoreload 2
%load_ext line_profiler

from masstodon.read.npy  import spectrum_from_npy
from masstodon.masstodon import masstodon_single, load_masstodon, masstodon_batch

data_path         = '/Users/matteo/Projects/masstodon/data/PXD001845/numpy_files/20141202_AMB_pBora_PLK_10x_40MeOH_1FA_OT_120k_10uscans_928_ETciD_8ms_15SA_19precZ/1'
mz, intensity     = spectrum_from_npy(data_path)
fasta             = "GAASMMGDVKESKMQITPETPGRIPVLNPFESPSDYSNLHEQTLASPSVFKSTKLPTPGKFRWSIDQLAVINPVEIDPEDIHRQALYLSHSRIDKDVEDKRQKAIEEFFTKDVIVPSPWTDHEGKQLSQCHSSKCTNINSDSPVGKKLTIHSEKSD"
charge            = 19
min_prob          = .8
isotopic_coverage = .999
std_cnt           = 3

m = masstodon_single(mz, intensity, fasta, charge,
                     orbitrap = True,
                     isotopic_coverage = isotopic_coverage,
                     min_prob = min_prob, 
                     std_cnt  = std_cnt,
                     deconvolution_method = 'nnls')

m.plotly()

m.plotly("/Users/matteo/Projects/masstodon/masstodon/dump/test_plotly", show=True)
m.imperator.dump_fig_to_json("/Users/matteo/Projects/masstodon/masstodon/dump/test_plotly/orbi.json")

