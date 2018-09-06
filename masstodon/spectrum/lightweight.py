from bisect import bisect
import numpy as np


class LightweightSpectrum(object):
    def __init__(self, min_mz, max_mz, total_intensities):
        self.clust_no = len(min_mz)
        self._spec = np.zeros(dtype = float,
                              shape = (self.clust_no*2,))
        # real spectrum
        self._spec[0::2] = min_mz
        self._spec[1::2] = max_mz
        self.intensity = total_intensities

    def __getitem__(self, key):
        i = bisect(self._spec, key)
        i_div, i_mod = divmod(i, 2)
        return i_div if i_mod else -1

    def __repr__(self):
        return self._spec.__repr__()

    def min_max_mz(self, i):
        return tuple(self._spec[(2*i):(2*i+2)])


def lightweight_spectrum(min_mz, max_mz, total_intensities):
    return LightweightSpectrum(min_mz,
                               max_mz,
                               total_intensities)


def test_lightweight_spectrum():
    L = lightweight_spectrum([0,5,8], 
                             [1,6,9],
                             [3,4,5])
    assert L[-2] == -1
    assert L[0]  ==  0
    assert L[.5] ==  0
    assert L[1]  == -1
    assert L[1.5]== -1
    assert L[4]  == -1
    assert L[5]  ==  1
    assert L[5.5]==  1
    assert L[6]  == -1
    assert L[10] == -1