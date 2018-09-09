import numpy as np


def by_rounding(x, y, digits):
    x_r    = np.around(x, digits)
    x_r, i = np.unique(x_r, return_inverse=True)
    y_agg  = np.bincount(i, weights=y)
    return x_r, y_agg


# def test_by_rounding():
#     x = np.array([.1111, .1112, .1131])
#     y = np.array([1,2,4])
#     X, Y = by_rounding(x, y, digits = 3)
#     assert all(X == np.array([.111, .113]))
#     assert all(Y == np.array([3,    4])