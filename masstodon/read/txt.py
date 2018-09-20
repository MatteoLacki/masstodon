import numpy as np


def spectrum_from_txt(path):
    """Read spectrum from a text file.

    Parameters
    ----------
    path : str
        Path to the mass spectrum file (mzxml, mzml, txt).
    Returns
    -------
    tuple : m/z ratios and intensities.
    """
    mz = []
    intensity = []
    with open(path) as f:
        for line in f:
            line = line.split()
            mz.append(float(line[0]))
            intensity.append(float(line[1]))
    return np.array(mz), np.array(intensity)