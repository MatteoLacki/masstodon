from collections import Counter, namedtuple
from math import floor, log10, inf
import numpy as np


from masstodon.models.polynomial import polynomial
from masstodon.models.spline import spline
from masstodon.spectrum.peak_clustering import \
    bitonic_clustering,\
    iter_cluster_ends,\
    min_diff_clustering
from masstodon.spectrum.lightweight import lightweight_spectrum
from masstodon.stats.gaussian import mean, sd, skewness
from masstodon.parse.threshold import parse as parse_threshold


class PeakClustering(object):
    """Common interface for peak clusterings."""
    def __init__(self):
        self.clusters = None

    def _iter_cluster_ends(self):
        yield from iter_cluster_ends(self.clusters)

    def __iter__(self):
        for s, e in self._iter_cluster_ends():
            yield self.x[s:e], self.w[s:e]

    def __len__(self):
        """Return the number of clusters."""
        return self.clusters[-1] + 1

    def left_ends_and_diffs(self):
        """Get left ends and widths of the clusters of basic peaks."""
        # the total number of diffs within clusters
        clusters_no = self.clusters[-1] - self.clusters[0]
        diffs_no    = len(self.x) - clusters_no - 1
        lefts       = np.zeros(shape=(diffs_no,), dtype=float)
        diffs       = lefts.copy()
        i_ = _i = 0
        for s, e in self._iter_cluster_ends():
            x = self.x[s:e]
            _i += len(x) - 1
            lefts[i_:_i] = x[:-1]
            diffs[i_:_i] = np.diff(x)
            i_ = _i
        return lefts, diffs


class MinDiff(PeakClustering):
    """The minimal difference peak clustering."""
    def fit(self, x, w, min_mz_diff = 1.1):
        self.x = x
        self.w = w
        self.clusters = min_diff_clustering(self.x,
                                            min_mz_diff)

def min_diff_clust(x, w, min_mz_diff = 1.1):
    mdc = MinDiff()
    mdc.fit(x, w, min_mz_diff)
    return mdc


# Problem: with absolute threshold, consecutive groups can overlay?
class Groups(object):
    def __init__(self):
        self.min_mz = []
        self.max_mz = []
        self.mean_mz = []
        self.sd_mz = []
        self.skewness_mz = []
        self.count = []
        self.intensity = []

    def get_stats(self,
                  local_spectra,
                  threshold,
                  out_trivial_intervals=True,
                  intensity_cumulant=sum):
        r_prev  = -inf
        thr, thr_type = parse_threshold(threshold)
        assert thr > 0.0, "Provide non-zero threshold."
        for local_mz, local_intensity in local_spectra:
            mean_mz = mean(local_mz, local_intensity)
            sd_mz   = sd(local_mz, local_intensity, mean_mz)
            if thr_type == "rel":
                l = mean_mz * (1.0 - thr)
                r = mean_mz * (1.0 + thr)
            else:
                l = mean_mz - thr
                r = mean_mz + thr
            # elif threshold_type == "precision":
            #     l = min(local_mz)
            #     r = max(local_mz)
            if l < r_prev: # overlaying intervals
                self.max_mz[-1] = l = (l+r_prev)/2.0
            if out_trivial_intervals and r > l:
                r_prev = r
                self.min_mz.append(l)
                self.max_mz.append(r)
                self.mean_mz.append(mean_mz)
                self.sd_mz.append(sd_mz)
                self.skewness_mz.append(skewness(local_mz, local_intensity, mean_mz, sd_mz))
                self.count.append(len(local_mz))
                self.intensity.append(intensity_cumulant(local_intensity))
        for k, v in self.__dict__.items():
            self.__dict__[k] = np.array(v)


class Bitonic(PeakClustering):
    """Clustering based on the bitonicity of intensities."""
    def __iter__(self):
        for s, e in self._iter_cluster_ends():
            yield self.x[s:e], self.w[s:e]

    def fit(self, x, w, 
                  min_mz_diff  = .15,
                  abs_perc_dev = .2):
        self.x = x
        self.w = w
        self.clusters = bitonic_clustering(self.x,
                                           self.w,
                                           min_mz_diff,
                                           abs_perc_dev)

    def get_stats(self,
                  threshold,
                  out_trivial_intervals=True):
        self.groups = Groups()
        self.groups.get_stats(self,
                              threshold,
                              out_trivial_intervals)
 
    def get_lightweight_spectrum(self):
        return lightweight_spectrum(self.groups.min_mz,
                                    self.groups.max_mz,
                                    self.groups.intensity)

    def get_smallest_diff_digits(self):
        """Return the number of digits of the smallest difference within the first bitonic cluster."""
        x, _ = next(self.__iter__())
        min_diff = np.diff(x)[0]
        return abs(floor(log10(min_diff)))

    def fit_diff_model(self, 
                       model=spline,
                      *model_args,
                     **model_kwds):
        lefts, diffs = self.left_ends_and_diffs()
        self.diff_model = model(lefts,
                                diffs,
                                *model_args,
                                **model_kwds)

    def fit_sd_model(self,
                     model = polynomial,
                     fit_to_most_frequent = True,
                    *model_args,
                   **model_kwds):
        # min_mz, max_mz, means, sds, skewnesses, counts, total_intensities, spreads = self.get_stats()
        if fit_to_most_frequent:
            cnts, freq    = list(zip(*Counter(self.groups.count).items()))
            self.sd_mz_c  = self.groups.count == cnts[np.argmax(freq)]
            self.sd_model = model(self.groups.mean_mz[self.sd_mz_c],
                                  self.groups.sd_mz[self.sd_mz_c])
        else:
            self.sd_model = model(self.groups.mean_mz, self.groups.sd_mz)

    # TODO: this should contain more class-specific code
    def plot_sd(self,
                plt_style = 'dark_background',
                show      = True):
        self.sd_model.plot(plt_style = plt_style,
                           show      = show)


def bitonic_clust(x,
                  w,
                  threshold,
                  min_mz_diff          = .15,
                  abs_perc_dev         = .2,
                  out_trivial_intervals= True,
                  model_diff           = spline,
                  model_diff_args      = [],
                  model_diff_kwds      = {},
                  fit_to_most_frequent = True,
                  model_sd             = polynomial,
                  model_sd_args        = [],
                  model_sd_kwds        = {}):
    bc = Bitonic()
    bc.fit(x, w, min_mz_diff, abs_perc_dev)
    bc.get_stats(threshold, out_trivial_intervals)
    bc.get_lightweight_spectrum()
    if model_diff:
        bc.fit_diff_model(model_diff,
                         *model_diff_args,
                        **model_diff_kwds)
    if model_sd:
        bc.fit_sd_model(model_sd,
                        fit_to_most_frequent,
                       *model_sd_args,
                      **model_sd_kwds)
    return bc

