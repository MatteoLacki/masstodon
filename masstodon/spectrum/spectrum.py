import numpy as np

from masstodon.read.npy           import spectrum_from_npy
from .orbitrap  import OrbitrapSpectrum
from .threshold import ThresholdSpectrum


def spectrum(mz        = np.array([]),
             intensity = np.array([]),
             min_mz_diff = 1.1,
             orbitrap  = False,
             threshold = 0.0,
             sort      = True):
    if orbitrap:
        spec = OrbitrapSpectrum(mz, intensity, sort)
        spec.min_mz_diff_clustering(min_mz_diff)
        spec.bitonic_clustering()
    else:
        spec = ThresholdSpectrum(mz, intensity, threshold, sort)
        spec.min_mz_diff_clustering(min_mz_diff)
    return spec


def load_spectrum(path):
    mz, intensity = spectrum_from_npy(path)
    spec = spectrum(mz,
                    intensity,
                    sort = False,
                    drop_duplicates = False)