"""Various operations on arrays."""
import numpy as np
import unittest


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
    if drop_duplicates:
        elements, counts = np.unique(x, return_counts = True)
        unique_entries   = np.isin(x, elements[counts == 1])
        x = x[unique_entries]
        y = y[unique_entries]
    if sort:
        i = np.argsort(x)
        x = x[i]
        y = y[i]
    return x, y


class Test_operations(unittest.TestCase):
    def test_dedup_sort(self):
        """Testing deduplication and sorting."""
        x = np.array([1, 2, 1, 3, 5])
        y = np.array(['a', 'b', 'c', 'd', 'e'])
        x_, y_ = dedup_sort(x, y)
        self.assertTrue(all(x_ == np.array([2,3,5])))
        self.assertTrue(all(y_ == np.array(['b','d','e'])))


if __name__ == "__main__":
    unittest.main()
