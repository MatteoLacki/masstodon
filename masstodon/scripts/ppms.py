%load_ext autoreload
%autoreload 2

from os.path import join as pjoin, exists as pexists
import matplotlib.pyplot as plt
import numpy as np


from masstodon.data.constants import infinity
from masstodon.masstodon import masstodon_single
from masstodon.read.mzml import read_mzml
from masstodon.plot.spectrum import plot_spectrum
from masstodon.spectrum.orbitrap import OrbitrapSpectrum

common_path = "/Users/matteo/Projects/masstodon/data/Belgian/2015_07/UBI_ORBI/12plus precursor"

min_prob          = .8
z                 = 12
isotopic_coverage = .999
std_cnt           = 3
timings           = True
ubiquitin         = 'MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG'

f = pjoin(common_path, "FRL_220715_ubi_714_ETD_0,03.mzXML")
x = next(read_mzml(f))
scan = int(x['num'])
mz = x['m/z array']
intensity = x['intensity array']

M, t = masstodon_single(mz, intensity, ubiquitin, z, "ubiquitin",
                        get_timings=True, threshold=5, threshold_type="ppm", orbitrap=True)
M.plotly("/Users/matteo/Projects/masstodon/spectrum.html")
list(M.ome.iter_molecule_estimates())


spec = OrbitrapSpectrum(mz, intensity, threshold=5, threshold_type="ppm")
# spec.plot()
spec.bitonic_clustering(groups_kwds={'ppm':5})
spec.bitonic_clustering()
# spec.plot_mz_diffs()
spec.bc.fit_diff_model()
# spec.bc.plot_sd()
spec.bc.groups.mean_mz

y = spec.bc.groups.max_mz - spec.bc.groups.min_mz
x = spec.bc.groups.mean_mz
plt.scatter(x,y)
plt.show()

np.percentile(spec.bc.groups.skewness_mz, np.linspace(0,100,10))
plt.hist(spec.bc.groups.sd_mz, 100)
plt.scatter(spec.bc.groups.mean_mz, spec.bc.groups.skewness_mz)
plt.show()
#TODO 
#   we have the mean. We need to simply modify the left and right ends, 
#   and that's it
