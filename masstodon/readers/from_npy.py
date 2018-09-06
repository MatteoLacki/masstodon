import numpy as np
import os

def spectrum_from_npy(data_path,
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
    mz =  np.load(os.path.join(data_path, mz_name))
    intensity =  np.load(os.path.join(data_path, intensity_name))
    return mz, intensity

