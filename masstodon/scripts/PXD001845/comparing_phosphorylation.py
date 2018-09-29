%load_ext autoreload
%autoreload 2
%load_ext line_profiler

from time import time

from masstodon.masstodon import Masstodon, masstodon_batch, masstodon_single, load_masstodon
from masstodon.ome.ome   import Ome
from masstodon.read.npy import spectrum_from_npy
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
    prec_mod['name'] += "_phosphorylated_" + str(i+1)
    prec_mod['modifications'] = {i+1: {'C_carbo': ptms['phosphorylation'].copy()}}
    precursors.append(prec_mod)

m = masstodon_single(mz, intensity, fasta, charge, '', 
                    orbitrap = True,
                    isotopic_coverage = isotopic_coverage,
                    min_prob = min_prob, 
                    std_cnt  = std_cnt)
m.plotlygl("/Users/matteo/Projects/masstodon/masstodon/dump", shape='triangles')
m.imperator.plot_ccs(plt_style = "default")

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
        prec_mod['name'] += "_phosphorylated_" + str(i+1)
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
n = masstodon_load("dump/phospho_search_targetted")


len(n.ome.G.nodes)
len(n.ome.G.edges)

import networkx as nx
from networkx.algorithms.flow import min_cost_flow


ome = m.ome

DG = nx.DiGraph()
demand = 0
for o in ome.observables():
    i = int(o.intensity)
    DG.add_node(o, demand=-i)
    demand += i
    for s in ome.G[o]:
        DG.add_edge(o,s,cost=len(s.modifications))
DG.add_node("sink", demand=demand)
for s in ome.sources():
    DG.add_edge(s,"sink")
res = min_cost_flow(DG, weight='cost')
minimal_intensity = {}
for s in ome.sources():
    s.intensity = res[s]["sink"]

for s in ome.sources():
    if s.intensity > 0:
        if s.modifications:
            i,p = list(s.modifications)[0]
            mes = "AA: {}\tAA_no: {}\tPTM: {}\tI: {}".format(s.fasta[i], i+1, s.modifications[(i,p)], s.intensity)
            
        else:
            mes = "non-modified\t\t\t\tI: {}".format(s.intensity)
        print(mes)


for o in ome.observables():
    o.name

