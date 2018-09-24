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

m = masstodon_single(mz, intensity, fasta, charge, '',
                     orbitrap = True,
                     isotopic_coverage = isotopic_coverage,
                     min_prob = min_prob, 
                     std_cnt  = std_cnt,
                     include_zero_intensities = True)
# m.dump('dump')
m.plotlygl("dump")

# n = masstodon_load("dump")
# n.plotlygl("dump", shape='triangles')




# trying out WebGL
import plotly
import plotly.graph_objs as go
import numpy as np

x = n.imperator.clust.groups.mean_mz
h = n.imperator.clust.groups.intensity.astype(int)
l = n.imperator.clust.groups.min_mz
r = n.imperator.clust.groups.max_mz


X = []
Y = []
for a,b,c in zip(l, r, h):
    X.append(a); Y.append(0)
    X.append(a); Y.append(c)
    X.append(b); Y.append(c)
    X.append(b); Y.append(0)

X = np.array(X)
Y = np.array(Y)

spectrum_bars = go.Scattergl(x = X, y = Y, name = "Peak Group")

fitted_int = []
fitted_mzs = []
fitted_to_int = []
for s in n.imperator.solutions:
    fitted_mzs.extend(list(s.mean_mz))
    fitted_int.extend(list(s.model.fitted()))
    fitted_to_int.extend(list(s.model.Y))
fitted_int      = np.array(fitted_int).astype(int)
fitted_mzs      = np.array(fitted_mzs)
fitted_to_int   = np.array(fitted_to_int).astype(int)
text_annotation = np.array(["fit {fit_int:.0f}<br>obs {fit2int:.0f}".format(fit_int=fit_int,
                                                                            fit2int=fit2int)
               for fit_int, fit2int in zip(fitted_int, fitted_to_int)])

fitted_dots = go.Scattergl(x = fitted_mzs,
                                     y = fitted_int,
                                     text      = text_annotation, # maps to labels
                                     hoverinfo = "text",          # show only labels
                                     mode      = "markers",       # default='lines'
                                     marker    = {"color" : "orange"},
                                     name      = "Fitted")

data = [spectrum_bars, fitted_dots]
layout = go.Layout(
    title = "Observed versus Fitted Spectrum",
    font  =dict(
        color="white"
    ),
    yaxis = dict(
        title = "Intensity",
        color = "white"
    ),
    xaxis = dict(
        title = "mass/charge [Th]",
        color = "white"
    ),
    plot_bgcolor = "black",
    paper_bgcolor= "black"
)
fig = go.Figure(data=data, layout=layout)
plotly.offline.plot(fig,
                    filename  = "dump/spectrum.html",
                    auto_open = True)

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
