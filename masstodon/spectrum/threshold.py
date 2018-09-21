from   collections import   namedtuple
import numpy       as       np
from   math        import   inf, floor, log10

from    masstodon.spectrum.lightweight import lightweight_spectrum

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


class ThresholdSpectrum(object):
    def fit(self, x, w, threshold):
        self.x = x
        self.w = w
        self.threshold = threshold
        self.clusters = np.arange(len(x))

    def get_smallest_diff_digits(self):
        """Return the number of digits of the smallest difference within the first bitonic cluster."""
        min_diff = min(self.groups.max_mz - self.groups.min_mz)
        return abs(floor(log10(min_diff)))

    def get_stats(self):
        l = []; r = []
        for p, c, n in smart_iter(self.x):
            l.append(c-min((c-p)/2.0, self.threshold))
            r.append(c+min((n-p)/2.0, self.threshold))
        l = np.array(l); r = np.array(r)
        self.groups = SimpleGroups(l, r, self.x, self.w)

    def get_lightweight_spectrum(self):
        self.ls = lightweight_spectrum(self.groups.min_mz,
                                       self.groups.max_mz,
                                       self.groups.intensity)


def threshold_spec(x, w, threshold):
    ts = ThresholdSpectrum()
    ts.fit(x, w, threshold)
    ts.get_stats()
    ts.get_lightweight_spectrum()
    return ts
