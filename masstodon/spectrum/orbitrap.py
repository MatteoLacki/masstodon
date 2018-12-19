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
import numpy   as np
from   bisect  import   bisect_left, bisect_right

from masstodon.spectrum.base    import Spectrum
from masstodon.spectrum.cluster import bitonic_clust
from masstodon.plot.spectrum    import plot_spectrum

class OrbitrapSpectrum(Spectrum):
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
                 threshold       = "0.0 Da",
                 sort            = True,
                 drop_duplicates = True,
                 drop_zeros      = True,
                 bc              = None,
                 mdc             = None):
        """Initialize the Spectrum."""
        super().__init__(mz, intensity, sort, drop_duplicates, drop_zeros, mdc)
        # parameters for spectra spawned as subspectra: needeed for convenience mainly.
        self.bc = bc
        self.threshold = threshold

    def __getitem__(self, interval):
        """Similar to filter, but return a class."""
        id_s = bisect_left(self.mz, interval.start)
        id_e = bisect_right(self.mz, interval.stop)
        return self.__class__(self.mz[id_s:id_e], 
                              self.intensity[id_s:id_e],
                              False,
                              False,
                              False,
                              self.bc,
                              self.mdc)

    def plot(self, 
             plt_style = "dark_background",
             peak_color= "white",
             clusters  = "",
             show      = True):
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
                                self.threshold,
                                min_mz_diff,
                                abs_perc_dev,
                                **kwds)

    def iter_subspectra(self, clustering):
        for s, e in clustering._iter_cluster_ends():
            yield self.__class__(mz              = self.mz[s:e],
                                 intensity       = self.intensity[s:e],
                                 sort            = False,
                                 drop_duplicates = False,
                                 drop_zeros      = False,
                                 bc              = self.bc.clusters[s:e],
                                 mdc             = self.mdc.clusters[s:e])

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

    def get_smallest_diff_digits(self):
        return self.bc.get_smallest_diff_digits()

    def get_lightweight_spectrum(self):
        return self.bc.get_lightweight_spectrum()

    def get_groups(self):
        return self.bc.groups
