from   math  import sqrt
import numpy as np


def mad(x, return_median=False):
    """Compute median absolute deviation from median.

    Args:
        x (np.array): preferably floats.
        return_median (boolean): should we return also the median?

    Return:
        out (float, or tupple): median absolute deviation for "x", potentially together with its median.
    """
    median = np.median(x)
    if return_median:
        return np.median(np.abs(x - median)), median
    else:
        return np.median(np.abs(x - median))


def mad_denoising(x, std_cnt = 100):
    """Denoise the signal using median absolute deviation.

    Simple robust approximation to standard deviation.
    https://en.wikipedia.org/wiki/Median_absolute_deviation
    """
    scaling   = 1.4826 # the magical scaling
    _mad, med = mad(x, return_median=True)
    sd        = _mad * scaling
    is_signal = np.abs(x - med) <= sd * std_cnt
    top       = med + sd * std_cnt
    bottom    = med - sd * std_cnt
    return is_signal, top, bottom
