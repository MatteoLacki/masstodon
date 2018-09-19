%load_ext autoreload
%autoreload 2
%load_ext line_profiler

from time import time

from masstodon.masstodon import Masstodon, masstodon_batch, masstodon_single, masstodon_load
from masstodon.ome.ome   import Ome
from masstodon.readers.from_npy import spectrum_from_npy

data_path     = '/Users/matteo/Projects/masstodon/data/PXD001845/numpy_files/20141202_AMB_pBora_PLK_10x_40MeOH_1FA_OT_120k_10uscans_928_ETciD_8ms_15SA_19precZ/1'
mz, intensity = spectrum_from_npy(data_path)

fasta             = "GAASMMGDVKESKMQITPETPGRIPVLNPFESPSDYSNLHEQTLASPSVFKSTKLPTPGKFRWSIDQLAVINPVEIDPEDIHRQALYLSHSRIDKDVEDKRQKAIEEFFTKDVIVPSPWTDHEGKQLSQCHSSKCTNINSDSPVGKKLTIHSEKSD"
charge            = 24
min_prob          = .8
isotopic_coverage = .999
std_cnt           = 3

# precursors = [{'fasta': fasta,
#                'name' : 'Abudhabizine',
#                'q'    : 24},
#                {'fasta': fasta[0:10],
#                'name' : 'Munichoson',
#                'q'    : 10}]

precursors = [{'fasta': fasta,
               'name' : '',
               'q'    : 24}]

# m = masstodon_batch(mz, intensity, precursors, 
#                     isotopic_coverage = isotopic_coverage,
#                     min_prob = min_prob, 
#                     std_cnt  = std_cnt)
# m = masstodon_single(mz, intensity, fasta, charge, '', 
#                     isotopic_coverage = isotopic_coverage,
#                     min_prob = min_prob, 
#                     std_cnt  = std_cnt)
# m.dump('dump')

n = masstodon_load("dump")
n.imperator.plotly_solutions('spectrum.html')


# %lprun -f masstodon_batch masstodon_batch(mz, intensity, precursors, isotopic_coverage=isotopic_coverage, min_prob=min_prob, std_cnt=std_cnt)
# n = masstodon_load("dump")
# o = masstodon_single(mz, intensity, fasta, charge, "Yani")
# p = masstodon_single(mz, intensity, fasta, charge, "Yani",
#                      deconvolution_graph_path="dump/deconvolution_graph.gpickle")

# from masstodon.precursor.precursor import precursor
# import intervaltree as iTree
# m = Masstodon()
# m.set_spectrum(mz, intensity)
# m.set_isotopic_calculator()
# t0 = time()
# mols = frozenset(m for p_kwds in precursors 
#                    for m in precursor(**p_kwds).molecules())
# emp_tree = iTree.IntervalTree()
# for subspec in m.subspectra:
#     s, e = subspec.interval
#     emp_tree[s:e] = subspec
# good_mols = []
# good_subspectra = set([])
# for mol, name in mols:
#     s, e = mol.interval(std_cnt = std_cnt)
#     touched_spectra = emp_tree[s:e]
#     if touched_spectra:
#         status = 'ok'
#         good_mols.append(mol)
#         good_subspectra |= touched_spectra
# t1 = time()
from masstodon.estimates_matcher.cz_simple import SimpleCzMatch



# precursors = [
#     dict(fasta = "AAAACCCKKK",
#          name  = "a-zine",
#          q     = 1),
#     dict(fasta = "VTRSQLM",
#          q     = 2,
#          name  = "b-zine"),
#     dict(fasta = "QQQQNNNKKKHIFFDDD",
#          q     = 2,
#          name  = "c-zine"),
# ]
# molecules = [
#     dict(name   = 'saliva',
#          formula= 'C100H202',
#          q      =  1),
#     dict(name   = 'shit',
#          formula= 'C102H205',
#          q      =  2),
#     dict(name   = 'lipstick',
#          formula= 'C100H202',
#          q      =  2),
#     dict(name   = "startrek's leather underwear",
#          formula= 'C100H202',
#          q      =  2)
# ]
