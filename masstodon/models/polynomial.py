"""Class that wraps up the polynomial fitting."""

import numpy as np

from masstodon.models.two_dimensional import Model2D
from masstodon.stats.descriptive      import mad_denoising


def fit_polynomial(x, y, degree = 3, coefs = True):
    _coefs = np.polyfit(x, y, degree)
    polynomial = np.poly1d(_coefs)
    if coefs:
        return polynomial, _coefs
    else:
        return polynomial


class Polynomial(Model2D):
    def __init__(self, degree = 3):
        """Initialize the Polynomial class.

        Parameters
        ----------
        degree : int
            The degree of the fitted polynomial.
        """
        assert degree >= 1, "The degree must be greater or equal to 1."
        self.degree = int(degree)

    def fit(self, x, y):
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
        self.x, self.y = x, y
        self.x_min = min(self.x)
        self.x_max = max(self.x)
        self._polynomial, self._coefs = fit_polynomial(self.x,
                                                       self.y,
                                                       self.degree,
                                                       True)

    def denoise_refit(self,
                    denoiser,
                  **denoiser_args):
        self.is_signal, _, _ = denoiser(self.errors(), **denoiser_args)
        x_s = self.x[self.is_signal]
        y_s = self.y[self.is_signal]
        self._polynomial, self._coefs = fit_polynomial(x_s,
                                                       y_s,
                                                       self.degree,
                                                       True)

    def coef(self):
        """Get the coefficients of the polynomial fitted with least squares."""
        return self._coefs


    def __call__(self, x):
        """Extrapolate the value of the fitted polynomial at 'x'.

        Parameters
        ----------
        x : np.array
            Values at which to extrapolate the polynomial.

        """
        return self._polynomial(x)

    def __repr__(self):
        return "Ich bin ein Polynom. Hast du Angst?"


def polynomial(x, y,
               degree        = 3,
               denoise_refit = {'denoiser': mad_denoising,
                                'std_cnt' : 100}):
    """Fit a scipy-spline to the data.

    Parameters:
    x : np.array
        1D control.
    y : np.array
        1D response.
    degree : int
        The degree of the fitted polynomial.
    """
    p = Polynomial(degree)
    p.fit(x, y)
    if denoise_refit:
        p.denoise_refit(**denoise_refit)
    return p