%load_ext autoreload
%autoreload 2
%load_ext line_profiler

from itertools import islice
import json
import os
from os.path import join as pjoin, exists as pexists, split as psplit
from sys import platform
import multiprocessing as mp
import pandas as pd
import numpy as np

from masstodon.scripts.PXD001845.iter_folders import iter_scans
from masstodon.spectrum.orbitrap import OrbitrapSpectrum
from masstodon.spectrum.cluster import Bitonic, Groups
from masstodon.stats.gaussian import mean, sd, skewness
from masstodon.spectrum.peak_clustering import bitonic_clustering, bitonic_iterator, fix_local_clustering
from masstodon.models.polynomial import polynomial
from masstodon.stats.descriptive import mad_denoising


# data paths
class WrongSystem(Exception):
    pass

if platform == "darwin":
    # check if you have your latte on skimmed soya milk with you
    data_path = "/Users/matteo/Projects/masstodon/data/PXD001845/numpy_files/"
    dump_path = "/Users/matteo/Projects/masstodon/dumps/many_processes/"
    processes_no = 4
elif platform == "linux":
    # check if you have long dirty hair
    data_path = "/home/matteo/masstodon/review_answer/numpy_files/"
    dump_path = "/mnt/disk/masstodon/dumps/many_processes/"
    processes_no = 1
elif "win" in platform:
    # don't check anything. no use.
    data_path = "C:/"
    raise WrongSystem(":D")
else:
    raise WrongSystem("Path not specified correctly.")


def iter_models():
    for mz, intensity, q, path in iter_scans(data_path):
        spec = OrbitrapSpectrum(mz, intensity)
        bc = Bitonic()
        bc.x = spec.mz
        bc.w = spec.intensity
        spec.bc = bc
        bc.clusters = bitonic_clustering(bc.x, bc.w, .15, .2)
        bc.fit_diff_model(model = polynomial,
                          degree= 2,
                          denoise_refit = {'denoiser': mad_denoising,
                                           'std_cnt' : 100})
        coefs = bc.diff_model.coef()
        deg   = len(coefs)
        row = {"C_" + str(deg - i - 1): v for i, v in enumerate(coefs)}
        row[q] = q
        row["path"] = path
        RANGE = np.arange(0,101,10)
        deciles = np.percentile(np.abs(bc.diff_model.res()), RANGE)
        for dec, val in zip(RANGE, deciles):
            row["perc_{dec}".format(dec=dec)] = val
        yield row

models_stats = list(islice(iter_models(), 100))
# models_stats.to_csv(pjoin(dump_path, "models_stats"), index=False)
pd.DataFrame(models_stats).to_csv(pjoin(dump_path, "models_stats"), index=False)
