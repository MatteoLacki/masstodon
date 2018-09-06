from __future__ import absolute_import, division, print_function
from IsoSpecPy.IsoSpecPy import IsoSpec as isospec

import numpy as np


def isospec_numpy(el_cnt,
                  el_mass,
                  el_prob,
                  threshold):
    """A wrapper around the IsoSpec software that neglects the configurations."""
    iso = isospec(el_cnt, 
                  el_mass,
                  el_prob,
                  threshold)
    mass, logprobability, _ = iso.getConfsRaw()
    mass = np.array(list(mass))
    prob = np.exp(list(logprobability))
    return mass, prob