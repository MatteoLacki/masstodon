from os.path import join as pjoin
from sys import platform

from masstodon.masstodon import masstodon_single, Masstodon
from masstodon.read.npy  import spectrum_from_npy
from masstodon.scripts.PXD001845.get_fasta import read_fasta
from masstodon.scripts.PXD001845.hecks_custom_ptms import modify_fasta
from masstodon.scripts.PXD001845.csv2fasta import csv2fasta

# data_path = "/Users/matteo/Projects/masstodon/data/PXD001845/numpy_files/20141202_AMB_pBora_AurA_10x_40MeOH_1FA_OT_60k_10uscans_924_EThcD_6ms_8CE_19precZ/15/"
# data_path = "/Users/matteo/Projects/masstodon/data/PXD001845/numpy_files/20141202_AMB_Bora_10x_40MeOH_1FA_OT_120k_10uscans_795_ETciD_4ms_15SA_22precZ/1/"
# data_path = "/Users/matteo/Projects/masstodon/data/PXD001845/numpy_files/20141202_AMB_Bora_10x_40MeOH_1FA_OT_120k_10uscans_728_ETD_4ms_24precZ/1/"
# data_path = "/Users/matteo/Projects/masstodon/data/PXD001845/numpy_files/20141202_AMB_Bora_10x_40MeOH_1FA_OT_120k_10uscans_728_ETD_4ms_24precZ/1/"

# czczmiel data
# data_path = "/Users/matteo/duch/20141202_AMB_Bora_10x_40MeOH_1FA_OT_120k_10uscans_728_ETD_4ms_24precZ/1/"
# data_path = "/home/matteo/masstodon/review_answer/numpy_files/20141202_AMB_Bora_10x_40MeOH_1FA_OT_120k_10uscans_728_ETD_4ms_24precZ/1/"

if platform == "darwin":
    # check if you have your latte on skimmed soya milk with you
    data_path = "/Users/matteo/Projects/masstodon/data/PXD001845/numpy_files/"
    dump_path = "/Users/matteo/Projects/masstodon/dumps/tests/"
elif platform == "linux":
    # check if you have long dirty hair
    data_path = "/home/matteo/masstodon/review_answer/numpy_files/"
    dump_path = "/mnt/disk/masstodon/dumps/tests/"

full_data_path = pjoin(data_path, "20141202_AMB_Bora_10x_40MeOH_1FA_OT_120k_10uscans_920_EThcD_8ms_15CE_19precZ/4/")
# folder = "20141202_AMB_Bora_10x_40MeOH_1FA_OT_120k_10uscans_728_ETD_4ms_24precZ"
folder = "20141202_AMB_Bora_10x_40MeOH_1FA_OT_120k_10uscans_920_EThcD_8ms_15CE_19precZ"

mz, intensity = spectrum_from_npy(full_data_path)
# from masstodon.plot.spectrum import plot_spectrum
# plot_spectrum(mz, intensity)

min_prob          = .8
isotopic_coverage = .999
std_cnt           = 3
q                 = 19

matching_csv = [f for f in csv2fasta if folder in f]
fasta = csv2fasta[matching_csv[0]]
fasta, modifications = modify_fasta(fasta)


m = masstodon_single(mz, intensity, fasta, q, '',
                     modifications = modifications,
                     orbitrap = True,
                     isotopic_coverage = isotopic_coverage,
                     min_prob = min_prob, 
                     std_cnt  = std_cnt,
                     include_zero_intensities = False)
m.dump(dump_path)
m.plotlygl(dump_path, show=False)
# m.plotlygl(dump_path)



# saving on czczmiel:
m.dump('/mnt/disk/masstodon/dumps/htmls/test')
m.plotlygl("/mnt/disk/masstodon/dumps/htmls/test", show=False)


# and now, the same with the method used in batch calculations:
data_path = "/Users/matteo/Projects/masstodon/data/PXD001845/numpy_files/20141202_AMB_Bora_10x_40MeOH_1FA_OT_120k_10uscans_795_ETciD_4ms_15SA_22precZ/1/"
data_path = "/Users/matteo/Projects/masstodon/data/PXD001845/numpy_files/20141202_AMB_Bora_10x_40MeOH_1FA_OT_120k_10uscans_795_ETciD_4ms_15SA_22precZ/"
mz, intensity     = spectrum_from_npy(data_path)
fasta             = "GAASMMGDVKESKMQITPETPGRIPVLNPFESPSDYSNLHEQTLASPSVFKSTKLPTPGKFRWSIDQLAVINPVEIDPEDIHRQALYLSHSRIDKDVEDKRQKAIEEFFTKDVIVPSPWTDHEGKQLSQCHSSKCTNINSDSPVGKKLTIHSEKSD"
charge            = 22
min_prob          = .8
isotopic_coverage = .999
std_cnt           = 3

# the mean of empty slice appeared here!
n = masstodon_single(mz, intensity, fasta, charge, '',
                     orbitrap = True,
                     isotopic_coverage = isotopic_coverage,
                     min_prob = min_prob, 
                     std_cnt  = std_cnt,
                     include_zero_intensities = True)
n.dump('dump')
n.plotlygl("dump")

