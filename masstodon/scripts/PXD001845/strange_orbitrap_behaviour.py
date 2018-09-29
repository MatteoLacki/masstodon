%load_ext autoreload
%autoreload 2
%load_ext line_profiler


from masstodon.masstodon import masstodon_single

from masstodon.scripts.PXD001845.get_fasta import read_fasta
from masstodon.scripts.PXD001845.hecks_custom_ptms import modify_fasta
from masstodon.scripts.PXD001845.hecks_spectrum import get_hecks_spectrum

# common
all_data  = "/Users/matteo/Projects/masstodon/data/PXD001845/numpy_files"
main_dump = "/Users/matteo/Projects/masstodon/masstodon/dump/debug"
folder    = "20141202_AMB_Bora_10x_40MeOH_1FA_OT_120k_10uscans_971_ETD_4ms_18precZ"

def fit(scan_no):
	mz, intensity = get_hecks_spectrum(all_data, folder, scan_no)
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

	m_dump = pjoin(main_dump, str(1))
	return m, tm, m_dump


m_scan = 1
m, tm, m_dump = fit(m_scan)
m.plotlygl(m_dump)


n_scan = 2
n, tn, n_dump = fit(n_scan)
n.plotlygl(n_dump)


