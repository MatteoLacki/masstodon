from   collections import   namedtuple
import numpy       as       np
from   math        import   inf, floor, log10

from masstodon.spectrum.orbitrap    import Spectrum
from masstodon.spectrum.lightweight import lightweight_spectrum


SimpleGroups = namedtuple('SimpleGroups',
                          'min_mz max_mz mean_mz intensity')
SimpleGroups.__new__.__defaults__ = tuple([] for _ in range(7))


def smart_iter(x):
    X = iter(list(x))
    prev_ = -inf
    curr_ = next(X)
    next_ = next(X)
    yield prev_, curr_, next_
    for next_ in X:
        prev_, curr_ = curr_, next_
        yield prev_, curr_, next_
    prev_, curr_ = curr_, next_
    yield prev_, curr_, inf


class ThresholdSpectrum(Spectrum):
    def __init__(self,   mz              = np.array([]),
                         intensity       = np.array([]),
                         threshold       = 0.0,
                         sort            = True,
                         drop_duplicates = True,
                         mdc             = None):
        assert threshold > 0.0, "Provide non-zero threshold."
        self.threshold = threshold
        super().__init__(mz, intensity, sort, drop_duplicates, mdc)
        l = []
        r = []
        for p, c, n in smart_iter(self.mz):
            l.append(c-min((c-p)/2.0, self.threshold))
            r.append(c+min((n-p)/2.0, self.threshold))
        l = np.array(l)
        r = np.array(r)
        self.groups = SimpleGroups(l, r, self.mz, self.intensity)

    def get_smallest_diff_digits(self):
        """Return the number of digits of the smallest difference within the first bitonic cluster."""
        min_diff = min(self.groups.max_mz - self.groups.min_mz)
        return abs(floor(log10(min_diff)))

    def get_lightweight_spectrum(self):
        return lightweight_spectrum(self.groups.min_mz,
                                    self.groups.max_mz,
                                    self.groups.intensity)

    def get_groups(self):
        return self.groups

    def iter_min_mz_diff_subspectra(self):
        for s, e in self.mdc._iter_cluster_ends():
            yield self.__class__(mz              = self.mz[s:e],
                                 intensity       = self.intensity[s:e],
                                 threshold       = self.threshold,
                                 sort            = False,
                                 drop_duplicates = True,
                                 mdc             = self.mdc.clusters[s:e])
