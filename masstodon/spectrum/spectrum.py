from   bisect            import bisect_left, bisect_right
try:
    import matplotlib.pyplot as plt
except RuntimeError:
    pass
import numpy             as np

from masstodon.arrays.operations import dedup_sort
from masstodon.measure.measure   import Measure
from masstodon.plot.spectrum     import plot_spectrum
from masstodon.spectrum.cluster  import min_diff_clust, bitonic_clust
from masstodon.stats.gaussian    import mean, sd, skewness


class Spectrum(Measure):
    """Prepare experimental spectrum.

    Parameters
    ----------
    spectrum : string, tuple, list, or Spectrum
        The path can end up with two extension:
        * txt, for spectra saved in tab separated format.
        * mzXml, for spectra saved with mxXml format.
        The tuple or list consist of two numpy arrays,
        one with m over z ratios and the other with intensities.
    mz_digits : float or int
        The number of digits after which the floats get rounded.
        E.g. if set to 2, then number 3.141592 will be rounded to 3.14.
    min_intensity : float
        Experimental peaks with lower height will be trimmed.
    percent_top_peaks : float
        Percentage of the heighest peaks in the spectrum to be included.
    """
    def __init__(self,
                 mz              = np.array([]),
                 intensity       = np.array([]),
                 sort            = True,
                 drop_duplicates = True,
                 only_positive   = True,
                 bc              = None,
                 mdc             = None,
                 mz_diff_model   = None,
                 min_mz_diff_bc  = None,
                 min_mz_diff_mdc = None):
        """Initialize the Spectrum."""
        self._store_names = ('m/z', 'intensity')
        self.clusters     = None
        self.mz, self.intensity = dedup_sort(mz,
                                             intensity,
                                             drop_duplicates,
                                             sort)
        self.mz        = self.mz[self.intensity > 0]
        self.intensity = self.intensity[self.intensity > 0]
        # parameters for spectra spawned as subspectra: needeed for convenience mainly.
        self.bc  = bc
        self.mdc = mdc
        self.mz_diff_model   = mz_diff_model
        self.min_mz_diff_bc  = min_mz_diff_bc
        self.min_mz_diff_mdc = min_mz_diff_mdc

    @property
    def mz(self):
        """Get mass over charge ratios"""
        return self.atoms

    @mz.setter
    def mz(self, mz):
        """Set m/z ratios."""
        self.atoms = mz

    @property
    def intensity(self):
        """Get intensities."""
        return self.masses

    @intensity.setter
    def intensity(self, intensity):
        """Set intensities."""
        self.masses = intensity

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
                              sort = False)

    @property
    def min_mz(self):
        return min(self.mz)

    @property
    def max_mz(self):
        return max(self.mz)

    @property
    def interval(self):
        return min(self.mz), max(self.mz)
 
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

    def plot(self, 
             plt_style = 'dark_background',
             peak_color= 'white',
             show      = True,
             clusters  = None):
        if clusters:
            if clusters is 'bitonic':
                clusters = self.bc.clusters 
            elif clusters is 'min_mz_diff':
                clusters = self.mdc.clusters
            else:
                raise KeyError("Wrong clustering name: should be either 'bitonic' or 'min_mz_diff'.")
        plot_spectrum(mz        = self.mz,
                      intensity = self.intensity,
                      clusters  = clusters,
                      plt_style = plt_style,
                      peak_color= peak_color,
                      show      = show)

    def bitonic_clustering(self,
                           min_mz_diff  = .15,
                           abs_perc_dev = .2,
                           **kwds):
        self.bc = bitonic_clust(self.mz,
                                self.intensity,
                                min_mz_diff,
                                abs_perc_dev,
                                **kwds)
        self.min_mz_diff_bc = min_mz_diff

    def min_mz_diff_clustering(self,
                               min_mz_diff = 1.1):
        if self.min_mz_diff_bc:
            assert min_mz_diff >= self.min_mz_diff_bc, "This can break up clusters: reduce min_mz_diff to levels lower than used for bitonic clustering."
        self.mdc = min_diff_clust(self.mz,
                                  self.intensity,
                                  min_mz_diff)
        self.min_mz_diff_mdc = min_mz_diff

    def iter_subspectra(self, clustering):
        for s, e in clustering._iter_cluster_ends():
            yield self.__class__(mz              = self.mz[s:e],
                                 intensity       = self.intensity[s:e],
                                 sort            = False,
                                 drop_duplicates = True,
                                 bc              = self.bc.clusters[s:e],
                                 mdc             = self.mdc.clusters[s:e],
                                 mz_diff_model   = self.mz_diff_model,
                                 min_mz_diff_bc  = self.min_mz_diff_bc,
                                 min_mz_diff_mdc = self.min_mz_diff_mdc)

    def iter_bitonic_subspectra(self):
        yield from self.iter_subspectra(self.bc)

    def iter_min_mz_diff_subspectra(self):
        yield from self.iter_subspectra(self.mdc)

    def plot_mz_diffs(self,
                      knots_no      = 1000,
                      plt_style     = 'dark_background',
                      all_diffs     = True,
                      cluster_diffs = True,
                      trend         = True,
                      show          = True):
        """Plot the clustering on the mass spectrum.

        Parameters
        ----------
            plt_style : str
                The type of the matplotlib style used in the plot.
            all_diffs : logical
                Plot all the differences, or only the ones in clusters?
            cluster_diffs : logical
                Plot the differences that belong to the cluster.
            trend : logical
                Plot the fitted trendline.
            show : logical
                Immediately show the plot? Alternatively, just add it to the present canvas.

        """
        plt.style.use(plt_style)
        if all_diffs:
            # Δ(m/z) = f(m/z)
            plt.scatter(self.mz[:-1],
                        np.diff(self.mz),
                        c = 'blue',
                        s = .5)
        if trend:
            self.bc.diff_model.plot(scatter_color = 'papayawhip',
                                    show = False)
        if show:
            plt.show()


def spectrum(mz        = np.array([]),
             intensity = np.array([]),
             sort      = True):
    spec = Spectrum(mz, intensity, sort)
    return spec
