%load_ext autoreload
%autoreload 2
%load_ext line_profiler

from time import time

from masstodon.masstodon import Masstodon, masstodon_batch, masstodon_single, masstodon_load
from masstodon.ome.ome   import Ome
from masstodon.readers.from_npy import spectrum_from_npy
from masstodon.data.ptms import ptms

dp = "/Users/matteo/Projects/masstodon/data/PXD001845/numpy_files/"
# data_path = dp+"20141202_AMB_Bora_10x_40MeOH_1FA_OT_120k_10uscans_920_ETD_4ms_19precZ/6"
data_path = dp+"20141202_AMB_pBora_PLK_10x_40MeOH_1FA_OT_120k_10uscans_928_ETciD_8ms_15SA_19precZ/1"
mz, intensity = spectrum_from_npy(data_path)

fasta             = "GAASMMGDVKESKMQITPETPGRIPVLNPFESPSDYSNLHEQTLASPSVFKSTKLPTPGKFRWSIDQLAVINPVEIDPEDIHRQALYLSHSRIDKDVEDKRQKAIEEFFTKDVIVPSPWTDHEGKQLSQCHSSKCTNINSDSPVGKKLTIHSEKSD"
charge            = 19
min_prob          = .8
isotopic_coverage = .999
std_cnt           = 3

# all phosphorylation sites
prec = {'fasta': fasta, 'q': charge, 'name': 'pBora'}
precursors = [prec]
for i, aa in enumerate(fasta):
    prec_mod = prec.copy()
    prec_mod['name'] += f"_phosphorylated_{i+1}"
    prec_mod['modifications'] = {i+1: {'C_carbo': ptms['phosphorylation'].copy()}}
    precursors.append(prec_mod)

# m = masstodon_single(mz, intensity, fasta, charge, '', 
#                     isotopic_coverage = isotopic_coverage,
#                     min_prob = min_prob, 
#                     std_cnt  = std_cnt)
m = masstodon_batch(mz, intensity, precursors, 
                    isotopic_coverage = isotopic_coverage,
                    min_prob = min_prob, 
                    std_cnt  = std_cnt)

m.dump('dump/phospho_search')
m.plotlygl("dump/phospho_search", shape='triangles')

# typical phosphorylation sites
prec = {'fasta': fasta, 'q': charge, 'name': 'pBora'}
precursors = [prec]
for i, aa in enumerate(fasta):
    if aa in ("S", "T", "Y"):
        prec_mod = prec.copy()
        prec_mod['name'] += f"_phosphorylated_{i+1}"
        prec_mod['modifications'] = {i+1: {'C_carbo': ptms['phosphorylation'].copy()}}
        precursors.append(prec_mod)
n = masstodon_batch(mz, intensity, precursors, 
                    isotopic_coverage = isotopic_coverage,
                    min_prob = min_prob, 
                    std_cnt  = std_cnt)

n.dump('dump/phospho_search_targetted')
n.plotlygl("dump/phospho_search_targetted", shape='triangles')

# 20141202_AMB_pBora_PLK_10x_40MeOH_1FA_OT_120k_10uscans_928_EThcD_6ms_10CE_19precZ/1
# strange peak plotting around 925.4-925.5
