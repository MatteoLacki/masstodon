"""TODO: document properly."""
try:
    import matplotlib.pyplot as plt
except RuntimeError:
    pass
import numpy as np

from masstodon.arrays.operations                 import dedup_sort
from masstodon.models.spline                     import spline
from masstodon.plot.spectrum                     import plot_spectrum
from masstodon.spectrum.orbitrap.peak_clustering import bitonic_clustering as bc,\
                                                        iter_cluster_ends

class Clustering(object):
    def __init__(self):
        self.clusters = None

    def fit(self, x, w):
        """Prepare self.clusters."""
        self.x = x
        self.w = w

    def __len__(self):
        """Get the number of clusters.

        Returns
        -------
            out : int
            The number of clusters.
        """
        return self.clusters[-1] + 1

    def iter_cluster_ends(self):
        """Iterate over consecutive cluster ends."""
        for s, e in iter_cluster_ends(np.nditer(self.clusters)):
            yield s, e

    def iter_clusters(self):
        for s, e in iter_cluster_ends(np.nditer(self.clusters)):
            yield self.x[s:e], self.w[s:e]

    def iter_x_and_x_diffs(self):
        """Iterate over m/z and m/z diffs.

        Iterate over tuples of arrays consisting of
        the left ends of spaces between consecutive
        peaks in a sequence of clusters.

        Yields:
        mz, mz_diff: tuples of np.arrays with mz values of all but the largest m/z in a cluster and the values of spaces between consecutive peaks.
        """
        for s, e in self.iter_cluster_ends():
            x      = self.x[s:e]
            x_diff = np.diff(x)
            yield x[:-1], x_diff

    def x_and_x_diffs(self):
        """Get x and x diffs."""
        diffs_no = len(self.x) - len(self)
        x = np.zeros(shape = (diffs_no,),
                     dtype = float)
        x_diffs = x.copy()
        i_ = _i = 0
        for x_l, x_diff_l in self.iter_x_and_x_diffs():
            _i += len(x_l)
            x[i_:_i]       = x_l
            x_diffs[i_:_i] = x_diff_l
            i_ = _i
        return x, x_diffs

    def fit_mz_diffs(self,
                     model=spline, 
                    *model_args,
                   **model_kwds):
        """Fit a spline to (m/z, Δm/z)."""
        x_lefts, x_diffs = self.x_and_x_diffs()
        self.x_diff_model = model(x_lefts, x_diffs, *model_args, **model_kwds)

    def x_diff(self, mz):
        """Return the fitted mz_diff for a given values of mz.

        Parameters
        ----------
        mz : np.array
            m/z values for which you need the values of the estimated m/z differences.
        """
        return self.x_diff_model(x)

    def plot(self,
             plt_style = 'dark_background',
             colors_no = 10,
             show      = True):
        """Plot the clustering on the mass spectrum.

        Parameters
        ----------
            plt_style : str
                The type of the matplotlib style used in the plot.
            colors_no : int
                The number of colors to cycle through.
            show : logical
                Immediately show the plot? Alternatively, just add it to the present canvas.

        """
        plot_spectrum(self.x,
                      self.w,
                      self.clusters,
                      plt_style,
                      colors_no,
                      show)


    def plot_mz_diffs(self,
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
            # plot all the subsequent m/z differences as function of m/z
            plt.scatter(self.x[:-1],
                        np.diff(self.x),
                        c = 'blue',
                        s = .5)
        if trend:
            self.x_diff_model.plot(plot_data = False,
                                   show      = False)
        if cluster_diffs:
            # mask those that are within clusters by bigger red dots.
            x_lefts, x_diffs = self.x_and_x_diffs()
            try:
                signal = self.x_diff_model.is_signal
                c = np.array(['papayawhip' if s else 'yellow' for s in signal])
            except AttributeError:
                c = 'papayawhip'
            plt.scatter(x_lefts,
                        x_diffs,
                        c = c,
                        s = 1.5)
        if show and (all_diffs or cluster_diffs):
            plt.show()


class BitonicClustering(Clustering):
    """Clustering based on bitonic intensities."""
    def __init__(self,
                 min_x_diff  = 0.15,
                 abs_perc_dev = 0.2):
        self.min_x_diff   = min_x_diff
        self.abs_perc_dev = abs_perc_dev

    def fit(self, x, w):
        self.clusters = bc(x,
                           w,
                           self.min_mz_diff,
                           self.abs_perc_dev)


def bitonic_clustering(x, w,
                       min_x_diff   = .15,
                       abs_perc_dev = .20,
                       diff_model   = spline, 
                      *diff_model_args, 
                     **diff_model_kwds):
    bc = BitonicClustering(min_x_diff, abs_perc_dev)
    bc.cluster(x, w, drop_duplicates = True, sort = True)
    bc.fit_mz_diffs(model = diff_model, *diff_model_args, **diff_model_kwds)
    return bc