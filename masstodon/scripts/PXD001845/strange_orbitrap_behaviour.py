%load_ext autoreload
%autoreload 2
%load_ext line_profiler


from masstodon.masstodon import masstodon_single

from masstodon.scripts.PXD001845.get_fasta import read_fasta
from masstodon.scripts.PXD001845.hecks_custom_ptms import modify_fasta
from masstodon.scripts.PXD001845.hecks_spectrum import get_hecks_spectrum

all_data = "/Users/matteo/Projects/masstodon/data/PXD001845/numpy_files"
folder   = "20141202_AMB_Bora_10x_40MeOH_1FA_OT_120k_10uscans_971_ETD_4ms_18precZ"
scan     = 1
mz, intensity = get_hecks_spectrum(all_data, folder, scan)

original_fasta 		 = read_fasta(folder)
fasta, modifications = modify_fasta(original_fasta)

charge            = 18
min_prob          = .8
isotopic_coverage = .999
std_cnt           = 3

m, tm = masstodon_single(mz, intensity, fasta, charge,
	                     orbitrap 		  	= True,
	                     isotopic_coverage 	= isotopic_coverage,
	                     min_prob 			= min_prob, 
	                     std_cnt  			= std_cnt,
	                     get_timings       	= True)

m.plotlygl("/Users/matteo/Projects/masstodon/masstodon/dump")
m.dump("dump")
