import numpy as np
import os
from os.path import join as pjoin, split as psplit, exists as pexists

path = "/Users/matteo/Projects/masstodon/data/PXD001845/numpy_files/"


def get_charge(experiment):
    return int(experiment.split("_")[-1].replace("precZ", ""))


def iter_scans(path,
               mz_file        = "mz.npy",
               intensity_file = "in.npy",
               results_folder = "res",
               paths_only     = False,
               stack_spectrum = False):
    """Iterate through the spectra saved to numpy objects and create folders for results.

    The presumed data storage folder scheme is:
        path/experiment/scan_no/mz_file
        path/experiment/scan_no/intensity_file
        path/experiment/scan_no/results_folder
    """
    for scan_dir, dirs, files in os.walk(path):
        try:
            mz_path           = pjoin(scan_dir, mz_file)
            intensity_path    = pjoin(scan_dir, intensity_file)
            mz                = np.load(mz_path)
            filepath, scan_no = psplit(scan_dir)
            experiment = psplit(filepath)[1]
            scan_no = int(scan_no)
            charge = get_charge(experiment)
            if not paths_only:
                intensity = np.load(intensity_path)   
                if stack_spectrum:
                    yield np.stack((mz, intensity), -1), charge, scan_dir
                else:
                    yield mz, intensity, charge, scan_dir
            else:
                yield mz_path, intensity_path, charge, scan_dir
            res_path = pjoin(scan_dir, results_folder)
            if not pexists(res_path):
                os.makedirs(res_path)
        except FileNotFoundError:
            pass
        except ValueError:
            pass


def non_modified_scans(path):
    for mz, intensity, charge, experiment in iter_scans(path):
        exp = experiment.split("/")[-2]
        if "AMB_Bora" in exp and "precZ" in exp and "ETD" in exp:
            yield mz, intensity, charge, experiment


