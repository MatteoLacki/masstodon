import pkg_resources
import json
import numpy as np

from masstodon.precursor.precursor import Precursor
from masstodon.spectrum.spectrum   import Spectrum


class Dataset(object):
    def __init__(self, precursor, spectrum, instrument):
        self.precursor = precursor
        self.spectrum = spectrum
        self.instrument = instrument

    def __repr__(self):
        out = "---- Dataset ----\n{}\n".format(self.precursor.__repr__())
        out += "{}Instrument{}\n".format(self.spectrum.__repr__(),
                                         self.instrument.__repr__())
        out += "----------------"
        return out


# TODO: this has to be modified to meet the old spectra specification.
# def get_dataset(dataset_name):
#     """Retrieve examplary spectra informations.

#     Parameters
#     ==========
#     dataset_name: string
#         Warning
#         =======
#         Can take values: substanceP or ubiquitin.
#     json : boolean
#         Should we read in json-saved file instead of python dictionary?

#     Returns
#     =======
#     A dictionary.
#     """
#     if dataset_name is 'substanceP':
#         from masstodon.data.substanceP import substanceP as mol
#     else:
#         raise AttributeError("No data set of that name.")
#     spectrum = Spectrum(spectrum=mol['spectrum'])

#     modifications = {int(k): v for k, v in
#                      mol['modifications'].items()}

#     precursor = Precursor(name=mol['name'],
#                           fasta=mol['fasta'],
#                           charge=mol['Q'],
#                           modifications=modifications,
#                           fragmentation_type="cz")

#     instrument = {}
#     if dataset_name == 'substanceP':
#         instrument['name'] = 'synapt'
#         instrument['wave height'] = 0
#         instrument['wave velocity'] = 300

#     elif dataset_name == 'ubiquitin':
#         instrument['name'] = 'orbitrap'
#         instrument['acquisition time'] = '10 ms'

#     return Dataset(precursor=precursor,
#                    spectrum=spectrum,
#                    instrument=instrument)
