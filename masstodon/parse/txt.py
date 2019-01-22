import numpy as np


def read_txt(path):
    """Read spectrum from an ASCI file.

    Parameters
    ----------
    path : str
        Path to the mass spectrum file (mzxml, mzml, txt).
        Every row in the file should contain m/z ratio separated from the intensity by some white-space character, such as space or tab.
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