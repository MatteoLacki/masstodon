%load_ext autoreload
%autoreload 2
from masstodon.spectrum.spectrum import spectrum
from masstodon.plot.spectrum import plot_spectrum
from masstodon.masstodon import masstodon_batch
import numpy as np

# path = "/Users/matteo/Downloads/erik/PFOA_experimental.txt"
path = "/Users/matteo/Downloads/erik/MSSD_UV0V_Org.txt"

with open(path,'r') as tsv:
    x = [line.strip().split('\t') for line in tsv]
masses = []
intensities = []
for y in x:
    masses.append(float(y[0]))
    intensities.append(float(y[1].replace(",", ".")))
intensities = np.array(intensities)
mz = np.array(masses)
# plot_spectrum(mz, intensities)


molecules = [{'formula':"C8H1F15O2",
              'q':q,
              'name':"erikium"} for q in range(1, 20)]

m, t = masstodon_batch(mz,
                       intensities * 10**7,
                       molecules                = molecules,
                       min_mz_diff              = .3,
                       orbitrap                 = False,
                       threshold                = .5,
                       get_timings              = True,
                       include_zero_intensities = False)

m.plotly(path="/Users/matteo/Downloads/erik/MSSD_UV0V_Org.html")
m.plotly(path="/Users/matteo/Downloads/erik/PFOA_experimental.html")
m.write(path="/Users/matteo/Downloads/erik/PFOA_experimental_res")

from masstodon.spectrum.spectrum import spectrum

s = spectrum(mz, intensities, threshold=1e-9)
s.plot()

