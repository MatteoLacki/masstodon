from pyteomics import mzml  # >= 3.41
from pyteomics import mzxml  # >= 3.41
import os

def read_mzml(path):
    """Read mzXML spectra.

    Generate a sequence of spectra from the path.
    
    Parameters
    ----------
    path : str
        Path to the mass spectrum file (mzxml, mzml, txt).
    Returns
    -------
    tuple : m/z ratios and intensities.
    """
    _, ext = os.path.splitext(path)
    ext = ext.lower()[1:]
    reader = {'mzxml': mzxml, 'mzml': mzml}[ext]
    with reader.read(path) as info:
        for spectrum in info:
            yield spectrum

