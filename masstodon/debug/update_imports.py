from masstodon.parse.spectrum_file import parse
import pandas as pd
import matplotlib.pyplot as plt

from masstodon.spectrum.spectrum import spectrum
from masstodon.spectrum.threshold import ThresholdSpectrum
from masstodon.parse.threshold import parse
from masstodon.spectrum.base import Spectrum

p = '/home/matteo/Projects/MassTodon/data/caroline/CRR_UMICH-ADH-5uM-550V-ECD-80V.csv'
D = pd.read_csv(p, skiprows=1)
mz = D['X(MassToCharge)']
intensity = D['Y(Counts)']

plt.plot(mz, intensity)
plt.show()


parse('0.1Da')
Spectrum(mz.values, intensity.values)


