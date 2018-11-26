try:
    import matplotlib.pyplot as plt
except RuntimeError:
    pass
try:
    import plotly
    import plotly.graph_objs as go
    from masstodon.plot.plotly import get_black_layout
    plotly_available = True
except ImportError:
    plotly_available = False

from   bisect  import   bisect_left, bisect_right
import numpy   as       np
from   os.path import   join as pjoin
import json

from masstodon.arrays.operations import dedup_sort
from masstodon.data.constants    import eps
from masstodon.plot.spectrum     import plot_spectrum
from masstodon.spectrum.cluster  import min_diff_clust
from masstodon.stats.gaussian    import mean, sd, skewness
from masstodon.misc.strings      import repr_long_list
from masstodon.spectrum.lightweight import lightweight_spectrum


class Spectrum(object):
    def __init__(self,   mz              = np.array([]),
                         intensity       = np.array([]),
                         sort            = True,
                         drop_duplicates = True,
                         drop_zeros      = True,
                         mdc             = None):
        """Initialize the Spectrum."""
        self.drop_duplicates = drop_duplicates
        self.drop_zeros      = drop_zeros
        if self.drop_zeros:
            mz        = mz[intensity > 0]
            intensity = intensity[intensity > 0]
        self.mz, self.intensity = dedup_sort(mz,
                                             intensity,
                                             self.drop_duplicates,
                                             sort)
        self.mdc = mdc

    def dump(self, path):
        np.save(pjoin(path, 'mz'), self.mz)
        np.save(pjoin(path, 'in'), self.intensity)

    def __repr__(self):
        """Represent the measure."""
        o = "{0}:\n\t{1} = {2}\n\t{3} = {4}\n".format(
              self.__class__.__name__,
              "m/z",
              repr_long_list(self.mz),
              "intensity",
              repr_long_list(self.intensity))
        return o

    def filter(self, mz_min, mz_max):
        """Filter part of the spectrum between mz_min and mz_max.

        Returns:
        tuple: m/z's and intensities.
        """
        id_s = bisect_left(self.mz, mz_min)
        id_e = bisect_right(self.mz, mz_max)
        return self.mz[id_s:id_e], self.intensity[id_s:id_e]

    def __getitem__(self, interval):
        """Similar to filter, but return a class."""
        id_s = bisect_left(self.mz, interval.start)
        id_e = bisect_right(self.mz, interval.stop)
        return self.__class__(self.mz[id_s:id_e], 
                              self.intensity[id_s:id_e],
                              False,
                              False,
                              False,
                              self.mdc)

    @property
    def min_mz(self):
        return min(self.mz)

    @property
    def max_mz(self):
        return max(self.mz)

    @property
    def interval(self):
        m = min(self.mz)
        M = max(self.mz)
        if m < M:
            return m, M
        else:
            return m * .999999, m * 1.000001

    def mean_mz(self):
        """Mean m/z weighted by intensities."""
        return mean(self.mz, self.intensity)

    def sd_mz(self):
        """Get standard deviation of m/z with probs induced by intensity."""
        return sd(self.mz, self.intensity)

    def skewness_mz(self):
        return skewness(self.mz, self.intensity)

    def total_intensity(self):
        return self.intensity.sum()

    def l1(self):
        """Get l1 norm."""
        return sum(self.intensity)

    def l2(self):
        """Get l2 norm."""
        return np.linalg.norm(self.intensity)

    def trim_intensity(self, cut_off):
        """Trim intensities below the provided cut off.

        Parameters
        ----------
        cut_off : float

        """
        self.trim(cut_off)

    def min_mz_diff_clustering(self,
                               min_mz_diff = 1.1):
        #TODO!!! Do something with that,
        # if self.min_mz_diff_bc:
        #     assert min_mz_diff >= self.min_mz_diff_bc, "This can break up clusters: reduce min_mz_diff to levels lower than used for bitonic clustering."
        self.mdc = min_diff_clust(self.mz,
                                  self.intensity,
                                  min_mz_diff)
        self.min_mz_diff_mdc = min_mz_diff

    def iter_min_mz_diff_subspectra(self):
        for s, e in self.mdc._iter_cluster_ends():
            yield self.__class__(mz              = self.mz[s:e],
                                 intensity       = self.intensity[s:e],
                                 sort            = False,
                                 drop_duplicates = False,
                                 drop_zeros      = False,
                                 mdc             = self.mdc.clusters[s:e])

    def get_min_mz_diff_subspectra(self):
        return list(self.iter_min_mz_diff_subspectra())

    def plot(self, 
             clusters  = None,
             plt_style = 'dark_background',
             peak_color= 'white',
             show      = True):
        if clusters:
            if clusters is 'min_mz_diff':
                clusters = self.mdc.clusters
            else:
                raise KeyError("Wrong clustering name: should be either 'bitonic' or 'min_mz_diff'.")
        plot_spectrum(mz        = self.mz,
                      intensity = self.intensity,
                      clusters  = clusters,
                      plt_style = plt_style,
                      peak_color= peak_color,
                      show      = show)

    def get_ligtweight_spectrum(self):
        raise NotImplementedError

    def get_smallest_diff_digits(self):
        raise NotImplementedError