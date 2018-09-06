# TODO: the name of the file does not match the nature of these functions that are 
# meaningful for any distribution.
import numpy as np
from math import sqrt, exp, log

from masstodon.data.constants import infinity

def max_intensity_mean(mz_g, intensity_g):
    """An idiotically simple estimator: the mass of the most intense peak."""
    return mz_g[np.argmax(intensity_g)]


def mean(mz_g, intensity_g):
    """Idiotically simple mean m/z estimator.

    Simply: a weighted mean of m/z."""
    return np.dot(mz_g, intensity_g)/sum(intensity_g)


def sd(mz_g, intensity_g, mean_mz = None):
    """Idiotically simple standard deviation of m/z estimator.

    Simply: a weighted sum of deviations to the mean.."""
    if mean_mz is None:
        mean_mz = mean(mz_g, intensity_g)
    probs = intensity_g/sum(intensity_g)
    return sqrt(np.dot((mz_g - mean_mz)**2, probs) )


def skewness(mz_g, intensity_g, mean_mz = None, sd_mz = None):
    """Idiotically simple estimator of the skewness of m/z for a group of peaks."""
    if len(intensity_g) >= 2:
        if mean_mz is None:
            mean_mz = mean(mz_g, intensity_g)
        if sd_mz is None:
            sd_mz = sd(mz_g, intensity_g, mean_mz)
        if sd_mz > 0:
            probs = intensity_g/sum(intensity_g)
            # return np.dot(((mz_g - mean_mz))**3, probs) / sd_mz**3
            try:
                thirds = np.dot(((mz_g - mean_mz))**3, probs)
                return exp(-3*log(sd_mz) + log(abs(thirds))) * np.sign(thirds)
            except ValueError:
                print(thirds, sd_mz)
                return 0.0
        else:
            return infinity
    else:
        return infinity