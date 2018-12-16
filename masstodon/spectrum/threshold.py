from   bisect  import   bisect_left, bisect_right
from   collections import   namedtuple
import numpy       as       np
from   math        import   inf, floor, log10

from masstodon.spectrum.orbitrap    import Spectrum
from masstodon.spectrum.lightweight import lightweight_spectrum


SimpleGroups = namedtuple('SimpleGroups',
                          'min_mz max_mz mean_mz intensity')
SimpleGroups.__new__.__defaults__ = tuple([] for _ in range(7))


class ThresholdSpectrum(Spectrum):
    def __init__(self,   mz              = np.array([]),
                         intensity       = np.array([]),
                         threshold       = 0.0,
                         threshold_type  = "mmu",
                         sort            = True,
                         drop_duplicates = True,
                         drop_zeros      = True,
                         mdc             = None):
        assert threshold > 0.0, "Provide non-zero threshold."
        assert threshold_type in ("mmu", "ppm", "Da")
        self.threshold = threshold
        self.threshold_type = threshold_type
        if threshold_type == "mmu":
            thr_type = "absolute"
            threshold /= 1000.0
        elif threshold_type == "ppm":
            thr_type = "relative"
        else:
            thr_type = "absolute"
        super().__init__(mz, intensity, sort, drop_duplicates, drop_zeros, mdc)
        X = iter(list(self.mz))
        x = next(X)
        thr = threshold if thr_type == "absolute" else threshold*x*1e-6
        l = [x - thr]
        r = []
        for x_next in X:
            d = min(thr, (x_next - x)/2.0)
            r.append(x + d)
            l.append(x_next - d)
            x = x_next
            thr = threshold if thr_type == "absolute" else threshold*x*1e-6
        r.append(x + thr)
        l = np.array(l)
        r = np.array(r)
        self.groups = SimpleGroups(l, r, self.mz, self.intensity)

    def get_smallest_diff_digits(self):
        """Return the number of digits of the smallest difference within the first bitonic cluster."""
        min_diff = min(self.groups.max_mz - self.groups.min_mz)
        return abs(floor(log10(min_diff)))

    #TODO: this should be a method of Spectrum.
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
                                 threshold_type  = self.threshold_type,
                                 sort            = False,
                                 drop_duplicates = True,
                                 drop_zeros      = False,
                                 mdc             = self.mdc.clusters[s:e])

    def __getitem__(self, interval):
        """Similar to filter, but return a class."""
        id_s = bisect_left(self.mz, interval.start)
        id_e = bisect_right(self.mz, interval.stop)
        return self.__class__(self.mz[id_s:id_e], 
                              self.intensity[id_s:id_e],
                              self.threshold,
                              self.threshold_type,
                              False,
                              False,
                              self.mdc)


