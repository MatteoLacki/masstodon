import numpy as np


def cdata2numpyarray(x):
    """Turn c-data into a numpy array.

    Parameters
    ----------
    x : cdata table
        A table of cdata from cffi.
    Returns
    -------
    res : array
        A numpy array of numbers.
    """
    res = np.empty(len(x))
    for i in range(len(x)):
        res[i] = x[i]
    return res


def get_mean_and_variance(X, weights):
    """Get mean and variance of X."""
    X = np.array(X)
    probs = np.array(weights)
    probs = probs / sum(probs)
    average = np.dot(X, probs)
    variance = np.dot((X - average) ** 2, probs)
    return average, variance


def check_charges(q, g):
    """Assert q and g are integers, q is positive, g nonnegative."""
    assert isinstance(q, int) and q > 0, "q must be a positive integer."
    assert isinstance(g, int) and g >= -1, "g must be greater than -1."
