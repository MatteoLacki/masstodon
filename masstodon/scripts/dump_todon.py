%load_ext autoreload
%autoreload 2
%load_ext line_profiler


from masstodon.masstodon        import masstodon_base
from masstodon.readers.from_npy import spectrum_from_npy

data_path     = '/Users/matteo/Projects/masstodon/data/PXD001845/numpy_files/20141202_AMB_pBora_PLK_10x_40MeOH_1FA_OT_120k_10uscans_928_ETciD_8ms_15SA_19precZ/1'
mz, intensity = spectrum_from_npy(data_path)

fasta  = "GAASMMGDVKESKMQITPETPGRIPVLNPFESPSDYSNLHEQTLASPSVFKSTKLPTPGKFRWSIDQLAVINPVEIDPEDIHRQALYLSHSRIDKDVEDKRQKAIEEFFTKDVIVPSPWTDHEGKQLSQCHSSKCTNINSDSPVGKKLTIHSEKSD"
charge = 24
min_prob = .8
isotopic_coverage = .999
std_cnt = 3

todon = masstodon_base(mz, intensity, fasta, charge,
                       std_cnt           = std_cnt,
                       isotopic_coverage = isotopic_coverage)

todon.imperator.plot_ccs()
todon.cz.plot(plt_style='seaborn')


# import jsonpickle

# frozen_todon = jsonpickle.encode(todon)

# with open('test.pickle', 'w') as f:
#     f.write(frozen_todon)

# wahaha = jsonpickle.decode(frozen_todon)

