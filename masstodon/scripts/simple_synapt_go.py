"""This needs revisiting. But the question is rather more about using
different error function altogether."""

from masstodon.masstodon import masstodon_single
from masstodon.data.substance_p import substance_p
from masstodon.data.constants import infinity


mz, intensity = substance_p['spectrum']

threshold         = 0.025
fasta             = 'RPKPQQFFGLM'
modifications     = {'11': {'C_carbo': {'H': 1, 'N': 1, 'O': -1}}}
# modifications = {}
q                 = 3
min_prob          = .8
isotopic_coverage = .999
std_cnt           = 3

m = masstodon_single(mz, intensity, fasta, q,
                     # min_mz_diff   = infinity,
                     modifications = modifications,
                     orbitrap      = False,
                     threshold     = threshold,
                     isotopic_coverage  = isotopic_coverage,
                     min_prob           = min_prob, 
                     std_cnt            = std_cnt,
                     include_zero_intensities = False)
# prec = next(m.ome.sources())
# for n in m.ome.observables():
#     print(m.ome.G[prec][n]["name"])
#     print(n.formula.tex_with_charges(ce=True))
#     print()

m.dump("")
m.plotlygl("", show=True)

sol = m.imperator.solutions[12]
[(i, sol.mean_mz) for i, sol in enumerate(m.imperator.solutions)]

sol.plot()
sol.model.fitted()



connected_component = sol.cc
total_intensities = m.imperator.groups.intensity
min_mz = m.imperator.groups.min_mz
max_mz = m.imperator.groups.max_mz
mean_mz = m.imperator.groups.mean_mz


from    networkx.linalg.attrmatrix  import attr_matrix
import  numpy                       as     np
try:
    import matplotlib.pyplot as plt
except RuntimeError:
    pass
from masstodon.models.nnls import nnls


cc          = connected_component
X, ordering = attr_matrix(cc, edge_attr='prob')
ordering    = np.array(ordering)
mol_columns = ordering < 0
peak_rows   = ordering >= 0
X           = X[:,mol_columns][peak_rows,:]
idx    = ordering[peak_rows]
total_intensities = total_intensities[idx]
mz_s   = min_mz[idx]
mz_e   = max_mz[idx]
mean_mz= mean_mz[idx]
Y           = total_intensities
model = nnls(X, Y)
model._coef

plt.subplot(1, 2, 1)
model.plot(show=False)


Y1 = np.concatenate((Y, np.zeros(X.shape[1])))
x1 = 1.0 - np.array(X.sum(axis=0)).flatten()
X1 = np.concatenate((X, np.diag(x1)))
model1 = nnls(X1, Y1)
model1._coef

plt.subplot(1, 2, 2)
model1.plot()



from masstodon.deconvolve.deconvolve import DeconvolutionProblem


dp = DeconvolutionProblem()
dp.fit(cc, )
