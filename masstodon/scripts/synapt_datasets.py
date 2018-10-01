%load_ext autoreload
%autoreload 2
%load_ext line_profiler

import pandas   as pd
from os         import listdir, makedirs
from os.path    import join as pjoin, exists as pexists

from masstodon.read.txt     import spectrum_from_txt
from masstodon.masstodon    import masstodon_single
from masstodon.data.constants import infinity

common_path = "/Users/matteo/Projects/masstodon/data/Belgian/2013_05_01/SUBP"
folders = ("wave_velocity_300", "wave_height_1.5")

def iter_data(path, folders):
    """Iterate over existing spectra in different folders with common root at the path."""
    for f in folders:
        pf = pjoin(path, f)
        for cff in listdir(pf):
            cffp = pjoin(pf, cff)
            try:
                yield cff, spectrum_from_txt(cffp)
            except Exception as e:
                pass


dump_folder = "/Users/matteo/Projects/masstodon/dumps/belgian/synapt"

threshold         = 0.025
fasta             = 'RPKPQQFFGLM'
modifications     = {'11': {'C_carbo': {'H': 1, 'N': 1, 'O': -1}}}
q                 = 3
min_prob          = .8
isotopic_coverage = .999
std_cnt           = 3
timings           = True
include_zero_intensities = False

bad = []

def res_with_figs():
    for i, (cff, (mz, intensity)) in enumerate(iter_data(common_path, folders)):
        try:
            m, t = masstodon_single(mz, intensity, fasta, q,
                                 min_mz_diff        = infinity,
                                 modifications      = modifications,
                                 orbitrap           = False,
                                 threshold          = threshold,
                                 isotopic_coverage  = isotopic_coverage,
                                 min_prob           = min_prob, 
                                 std_cnt            = std_cnt,
                                 get_timings        = timings,
                                 include_zero_intensities = include_zero_intensities)
            df = pjoin(dump_folder, "{i}_".format(i=i) + cff)
            if not pexists(df):
                makedirs(df)
            m.dump(df)
            m.write(df)
            m.plotlygl(df, shape='rectangles', show=False)
            print(f'Finished with {cff}')
        except AssertionError:
            bad.append((cff, mz, intensity))

res_with_figs()
bad


def get_WH_WV(f):
    f = f.replace(" ","")
    f = f.split(".")[0]
    WH, WV = f.split("-")[3:5]
    WH = WH[2:]
    if WH == '1,5':
        WH = 1.5
    else:
        WH = int(WH)/100
    WV = int(WV[2:])
    return WH, WV

def iter_outcomes():
    for i, (cff, (mz, intensity)) in enumerate(iter_data(common_path, folders)):
        WH, WV = get_WH_WV(cff)
        m, t = masstodon_single(mz, intensity, fasta, q,
                             min_mz_diff        = infinity,
                             modifications      = modifications,
                             orbitrap           = False,
                             threshold          = threshold,
                             isotopic_coverage  = isotopic_coverage,
                             min_prob           = min_prob, 
                             std_cnt            = std_cnt,
                             get_timings        = timings,
                             include_zero_intensities = include_zero_intensities)
        row = {"i":i, "exp": cff, "WH": WH, "WV": WV}
        row.update(m.imperator.errors())
        row.update({"t_"+str(n): T for n,T in t})
        for s in ('ETDorHTR', 'ETnoD_PTR_fragments', 'ETnoD_precursor', 'PTR_precursor'):
            row["cz_simple."+str(s)] = int(m.cz_simple.intensities[s])
            row["cz."+str(s)]        = int(m.cz.intensities[s])
        for i in range(11):
            row["cz_simple.prob." + str(i)] = m.cz_simple.probabilities["fragmentation_bond"].get(i, 0.0)
            row["cz.prob." + str(i)] = m.cz.probabilities["fragmentation_bond"].get(i, 0.0)
        yield row
        print('Finished with {cff}'.format(cff=cff))

results_for_plot = pd.DataFrame(iter_outcomes())
results_for_plot.to_csv(
    "/Users/matteo/Projects/masstodon/dumps/belgian/synapt_results.csv",
    index = False)







# getting all the results again...
def iter_outcomes_raw():
    for i, (cff, (mz, intensity)) in enumerate(iter_data(common_path, folders)):
        WH, WV = get_WH_WV(cff)
        m,t = masstodon_single(mz, intensity, fasta, q,
                             min_mz_diff        = infinity,
                             modifications      = modifications,
                             orbitrap           = False,
                             threshold          = threshold,
                             isotopic_coverage  = isotopic_coverage,
                             min_prob           = min_prob, 
                             std_cnt            = std_cnt,
                             get_timings        = timings,
                             include_zero_intensities = include_zero_intensities)
        yield m, t, cff, WH, WV, mz, intensity

D = list(iter_outcomes_raw())
m, t, cff, WH, WV, mz, intensity = [d for d in D if d[3] == 1.5 and d[4] == 10][0]
m.plotlygl("dump")

prec = m.good_mols[3]
f = prec.formula

m.spec.total_intensity()
prec.intensity


 m,t = masstodon_single(mz, intensity, fasta, q,
                         min_mz_diff        = infinity,
                         modifications      = modifications,
                         orbitrap           = False,
                         threshold          = threshold,
                         isotopic_coverage  = isotopic_coverage,
                         min_prob           = min_prob, 
                         std_cnt            = std_cnt,
                         get_timings        = True)

m.plotlygl("dump")

i = 21
m.imperator.solutions[i].plot()

S = m.imperator.solutions[i]
S.idx
S.mz_s
S.mz_e
S.plot()
S.X[:,0]

for mol in m.good_mols:
    print(mol, mol.mean_mz)


mol_columns = np.array([N < 0  for N in S.cc])
peak_rows   = np.array([N >= 0 for N in S.cc])
X, ordering = attr_matrix(S.cc, edge_attr='prob')
X           = X[:,mol_columns][peak_rows,:]
ordering    = np.array(ordering)
S.idx    = ordering[ordering >= 0]
S.total_intensities = total_intensities[S.idx]
S.mz_s   = min_mz[S.idx]
S.mz_e   = max_mz[S.idx]
S.mean_mz= mean_mz[S.idx]
Y           = S.total_intensities
if include_zero_intensities:
    Y = np.concatenate((Y, np.zeros(X.shape[1])))
    x = 1.0 - np.array(X.sum(axis=0)).flatten()
    X = np.concatenate((X, np.diag(x)))
S.X = X
S.Y = Y
S.model = nnls(X, Y)



# get the precursor problem: ~ 674 Th
for i, sol in enumerate(m.imperator.solutions):
    print(i, sol.mean_mz)
# how to retrieve which molecules are here?
# get the over-estimated precursor
prec_obs = [o for o in m.ome.observables() if o.intensity > 9*10**8][0]
prec = [s for s in m.ome.G[prec]][0]

prec_obs.isotopologues().probability

# there must be something wrong with attr_matrix function
m.imperator.G[prec_obs]
m.imperator.G.nodes

for e in S2.cc.edges(data=True):
    print(e)


S2 = m.imperator.solutions[1]
S2.plot_fancy()
S2.X
S2.Y
S2.plot()
S2.model.plot()
S2.model.coef()
list(S2.iter_estimates())

from networkx import draw
import matplotlib.pyplot as plt
draw(S2.cc)
plt.show()


total_intensities = m.groups.intensity
min_mz = m.groups.min_mz
max_mz = m.groups.max_mz
mean_mz = m.groups.mean_mz


# mol_columns = np.array([N < 0  for N in S2.cc])
# peak_rows   = np.array([N >= 0 for N in S2.cc])
include_zero_intensities = False


X, ordering = attr_matrix(S2.cc, edge_attr='prob')
ordering    = np.array(ordering)
mol_columns = ordering < 0
peak_rows   = ordering >= 0
X           = X[:,mol_columns][peak_rows,:]
S2.idx      = ordering[peak_rows]
S2.total_intensities = total_intensities[S2.idx]
S2.mz_s   = min_mz[S2.idx]
S2.mz_e   = max_mz[S2.idx]
S2.mean_mz= mean_mz[S2.idx]
Y          = S2.total_intensities
if include_zero_intensities:
    Y = np.concatenate((Y, np.zeros(X.shape[1])))
    x = 1.0 - np.array(X.sum(axis=0)).flatten()
    X = np.concatenate((X, np.diag(x)))
S2.X = X
S2.Y = Y
X.shape
Y.shape
S2.model = nnls(X, Y)
S2.model.plot()

S2.model._coef
#  write a deconvolution test.






S2.cc.nodes()
[e for e in S2.cc.edges(data=True)]

attr_matrix(S2.cc)
WW = attr_matrix(S2.cc, edge_attr='prob')[0]
np.round(WW,2)
WW[4:6,]


mol_columns = np.array([N < 0  for N in S2.cc])
peak_rows   = np.array([N >= 0 for N in S2.cc])




X = np.delete(S2.X, (5), 0)
Y = np.delete(S2.Y, (5), 0)

scipy.optimize.nnls(X,Y, maxiter=2)

scipy.optimize.nnls(X[:,0],Y, maxiter=2)
scipy.optimize.nnls(X[:,1],Y, maxiter=2)

model = nnls(X, Y)
model.coef()
model._res
model._coef




