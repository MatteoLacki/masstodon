try:
    import matplotlib.pyplot as plt
except RuntimeError:
    pass
import numpy             as np
from   scipy.interpolate import LSQUnivariateSpline

from masstodon.arrays.operations      import dedup_sort
from masstodon.models.two_dimensional import Model2D
from masstodon.stats.descriptive      import mad_denoising


class Spline(Model2D):
    def fit(self, x, y,
            drop_duplicates = True,
            sort            = True):
        """Fit a spline to 2D data.

        Parameters
        ----------
        x : np.array
            One-dimensional control variables.
        y : np.array
            One-dimensional response variable.
        drop_duplicates : logical
            Drop duplicate rows in 'xy' data-frame.
        sort : logical
            Sort the 'xy' data-frame by 'x'.
        Returns
        -------
        tuple of np.arrays: fixed 'x' and 'y'.
        """
        self.x, self.y = dedup_sort(x, y, drop_duplicates, sort)
        t = self.x_percentiles(len(self.x)//1000)
        self.x_min = t[0]
        self.x_max = t[-1]
        t = t[1:-1]
        self._spline = LSQUnivariateSpline(self.x, self.y, t)


    def denoise_refit(self,
                      denoiser,
                    **denoiser_args):
        self.is_signal, _, _ = denoiser(self.res(), **denoiser_args)
        x_s = self.x[self.is_signal]
        y_s = self.y[self.is_signal]
        t = self.x_percentiles(len(x_s)//1000)[1:-1]
        self._spline = LSQUnivariateSpline(x_s, y_s, t)

    def __call__(self, x):
        return self._spline(x)


    def __repr__(self):
        return "Ich bin ein Spline. Hast du Angst?"


def spline(x, y, 
           drop_duplicates = True,
           sort            = True,
           denoise_refit   = {'denoiser': mad_denoising,
                              'std_cnt' : 100}):
    """Fit a scipy-spline to the data.

    Arguments:
        x : np.array
        y : np.array
    """
    s = Spline()
    s.fit(x, y, drop_duplicates, sort)
    if denoise_refit:
        s.denoise_refit(**denoise_refit)
    return s