import numpy as np
import os

def spectrum_from_npy(data_path      = '.',
                      mz_name        = 'mz.npy',
                      intensity_name = 'in.npy'):
    """Import spectrum from numpy files.

    It is assumed that two files exist:
    one with m/z values, the other containing the intensities.
    They are assumed to be in the same folder.

    Parameters
    ----------
    data_path : str
        The path to the folder containing the files.
    mz_name : str
        The name of the file that contains m/z values.
    intensity_name : str
        The name of the file that contains the intensity values.
    """
    mz = np.load(os.path.join(data_path, mz_name))
    intensity = np.load(os.path.join(data_path, intensity_name))
    return mz, intensity


def read_npy(path):
    file = os.path.basename(path)
    assert file in ('mz.npy', 'in.npy'), "Supply a path to some mz.npy or in.npy file. The other file should be in the same folder."
    if file == 'mz.npy':
        d_p = path.replace('mz.npy', '')
    elif file == 'in.npy':
        mz_p = path.replace('in.npy', '')

