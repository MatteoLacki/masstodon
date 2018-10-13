%load_ext autoreload
%autoreload 2
%load_ext line_profiler

from masstodon.read.npy  import spectrum_from_npy
from masstodon.masstodon import masstodon_single, load_masstodon, masstodon_batch

data_path         = '/Users/matteo/Projects/masstodon/data/PXD001845/numpy_files/20141202_AMB_pBora_PLK_10x_40MeOH_1FA_OT_120k_10uscans_928_ETciD_8ms_15SA_19precZ/1'
mz, intensity     = spectrum_from_npy(data_path)
fasta             = "GAASMMGDVKESKMQITPETPGRIPVLNPFESPSDYSNLHEQTLASPSVFKSTKLPTPGKFRWSIDQLAVINPVEIDPEDIHRQALYLSHSRIDKDVEDKRQKAIEEFFTKDVIVPSPWTDHEGKQLSQCHSSKCTNINSDSPVGKKLTIHSEKSD"
charge            = 19
min_prob          = .8
isotopic_coverage = .999
std_cnt           = 3

m = masstodon_single(mz, intensity, fasta, charge,
                     orbitrap = True,
                     isotopic_coverage = isotopic_coverage,
                     min_prob = min_prob, 
                     std_cnt  = std_cnt)

m.plotlygl("/Users/matteo/Projects/masstodon/masstodon/dump")
m.dump("dump")

# n = load_masstodon("dump")
# path = "dump"

# from masstodon.deconvolve.divide_ed_impera import load_imperator
# load_imperator('dump/deconvolution_graph.gpickle')
p = m.imperator.solutions[100]
p.plot()



import numpy as np
from scipy.optimize import linprog
# def quant_deconv(X, Y):

X = p.X
Y = p.Y
q = .25
l1= .00

K = Y.shape[0]
D = X.shape[1]

c = np.concatenate( ( (1-q) * np.ones(shape=(K,1)),
                        q   * np.ones(shape=(K,1)),
                        l1  * np.ones(shape=(D,1))) )
c.shape = (c.shape[0],)

A_ub = np.vstack((
    np.hstack((-np.identity(K),
                np.zeros(shape=(K,K)),
                X )),
    np.hstack(( np.zeros(shape=(K,K)),
               -np.identity(K),
               -X )) ))

b_ub  = np.concatenate((Y,-Y))
model = linprog(c    = c,
                A_ub = A_ub,
                b_ub = b_ub)

model['x'][-D:]

help(np.dot)

np.abs(Y - X @ model['x'][-D:])


plt.vlines(x = x, ymin=[0], ymax=self.Y, colors=lines_col)



