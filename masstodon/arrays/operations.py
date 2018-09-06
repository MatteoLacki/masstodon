"""Various operations on arrays."""

import numpy as np
import pandas as pd


def is_sorted(x):
    return np.all(np.diff(x))


def dedup_sort(x, y, 
               drop_duplicates = True,
               sort            = True):
    """Remove dulicate x entries in x and for the corresponding y indices. 
    Sort by x.

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
    if drop_duplicates or sort:
        d = pd.DataFrame({'x':x, 'y':y})
        if drop_duplicates:
            d = d.drop_duplicates(subset='x', keep=False)
        if sort:
            d = d.sort_values(['x'])
        return d.x.values, d.y.values
    else:
        return x, y
