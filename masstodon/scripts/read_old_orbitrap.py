%load_ext autoreload
%autoreload 2
%load_ext line_profiler

from masstodon.read.mzml import mzml

datapath = "/Users/matteo/Projects/masstodon/data/Belgian/2015_07/UBI_ORBI/12plus precursor/FRL_220715_ubi_714_ETD_0,03.mzXML"
it = mzml(datapath)
info = next(it)
infos = list(it)
len(infos)


mz        = info['m/z array']
intensity = info['intensity array']
mz        = mz[intensity > 0]
intensity = intensity[intensity > 0]

mz_old = mz
intensity_old = intensity
mz
mz_old

m0 = set(mz)
m1 = set(mz_old)

len(m0 and m1)
len(m0)
len(m1)

# what is common in most mzs?


from masstodon.plot.spectrum import plot_spectrum


plot_spectrum(mz, intensity, show=False)
plot_spectrum(mz_old, -intensity_old)


from masstodon.spectrum.spectrum import spectrum


spec = spectrum(mz, intensity)
spec.bitonic_clustering()
spec.plot()


