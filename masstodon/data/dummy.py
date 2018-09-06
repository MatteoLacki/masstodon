"""Some simple an inoccuous dummy data and how to plot it."""

import numpy as np

mz_dummy        = np.array([1.1, 1.2, 1.3, 1.4, 2.2, 2.3, 2.4, 2.5, 2.6, 2.8, 2.9, 3.0])
intensity_dummy = np.array([0.1, 1.2, 2.3, 1.4, 0.2, 2.3, 4.4, 2.5, 1.6, 0.8, 2.9, 0.2])


def plot_dummy():
    """Plot the dummy data."""
    from MassTodonPy.plotters.spectrum import plot_spectrum
    from MassTodonPy.Spectra.peak_clustering import mz_bitonic
    clusters_dummy = mz_bitonic(mz_dummy, intensity_dummy)
    plot_spectrum(mz_dummy, intensity_dummy, list(clusters_dummy))
