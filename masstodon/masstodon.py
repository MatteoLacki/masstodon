# -*- coding: utf-8 -*-
#   Copyright (C) 2016 Mateusz Krzysztof Łącki and Michał Startek.
#
#   This file is part of MassTodon.
#
#   MassTodon is free software: you can redistribute it and/or modify
#   it under the terms of the GNU AFFERO GENERAL PUBLIC LICENSE
#   Version 3.
#
#   MassTodon is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
#   You should have received a copy of the GNU AFFERO GENERAL PUBLIC LICENSE
#   Version 3 along with MassTodon.  If not, see
#   <https://www.gnu.org/licenses/agpl-3.0.en.html>.


# Here is our little pet MassTodon!
# ................................................................................
# ............~MMM:...............................................................
# ..........:M....NM..............................................................
# ..........8M,....M,.............................................................
# ............M~...8O.............................................................
# ..........+.MM,..,M..............................?MMMMMN=.......................
# ..........:..:M:..M..........................DMMM~..  ....MMMM$.................
# ..........+...DN..M...............+MMMMMMMMMMM............... $MMM=.............
# ..........?...OM..M=............ZMD...... .......................ZMM?...........
# ..........:...?M..8O...........M,...................................MM,.........
# ..........+...7M..~M..........M......................................OMM........
# ..........+...IM...M.........MI....................................... MM.......
# ..........:...ZM...M.........M...M?.8D..................................MM......
# ..........:...DM...M:.......M:..M....,O.................................?M:.....
# ..............NM...~M......OM...M....8:..................................MM.....
# ...........$..ZM....MO.....M.....DMMM.....................................M,....
# ...........Z..,M.....MM..8M:..............................................MO....
# ..........:Z...MN......DMO................................................MM....
# ..........:Z...~M~........................................................7M....
# ..........?Z....$M7........................................................M....
# ..........+Z.....=MM...........MZ..........................................,M...
# ..........~Z.......MMMM=....IMMMM..~M....................................:M.M...
# ..........=Z......:...IMMMMMZ..MI..NMM....................................NM=...
# ..........+Z......MMZ.........OM...M,.M....................................M+...
# ..........:Z.....Z.M.,MMMIIDMMM....M..M.N~..................................M...
# ..........:Z.....Z..MM...........~MD..M.MM................................M,M7..
# ..........:Z.....Z....7MMMO:~8MMM7...ND.MN................................MM:+..
# ..........:Z.....Z..................?I..MD................................M.M...
# ..........=$.....?.....................MZN................................M+....
# ..........ZZ.....7....................MD.M........MMM..N,..M...M7M........MO....
# ..........Z$~....?.......................M........M..M.~D.~N..:,.M........MO....
# ..........ZZ,....$.... .......=..........M........M:.M..D..M..+=.M,.......MO....
# ........~.ZZ?....$....=.......=..........MM......?M,.,...........MM......8M=....
# .$I.,.IZ,I.+,~I.,....~~..7Z.Z.:...=.,$Z...~8MMMMMMI..............~MMMMMMMM8.....
# ................................................................................
# ................................................................................
# ................................................................................
# .7O......O............................ZZZZZZZ..............Z....................
# .7......I$...............................Z.................Z....................
# .7.Z....O7....,ZZZ.....OZZ......ZZO .....Z.....OZZ......OZOZ....,ZZO....ZOOZ+...
# .7..$..Z.7...7....Z...Z....=.......O.....Z....O....?...=...?...I....Z...Z....O..
# .7..O..~.7........Z.. Z........7.........Z....O....O...O...O...Z....Z...O....Z..
# .7...ZZ..7.....8Z.Z.....O.......$Z.......Z....O....O...O...O...Z....Z...O....Z..
# .7...~...7...Z....Z.......7........O.....Z....O....O...O...O...Z....Z...O....Z..
# .7.......7...O....O...O... I.......Z.....Z....Z... $...O...O...$....O...O....Z..
# .7.......7....OOOOO....ZOZI....:ZOZ......Z.....ZOZ7.....OZZ.....7ZOZ....Z....Z..
# ................................................................................
import json
import numpy  as      np
from   os.path import join as pjoin
from   time   import  time

from masstodon.deconvolve.divide_ed_impera import divide_ed_impera, Imperator
from masstodon.estimates_matcher.cz        import CzMatch
from masstodon.estimates_matcher.cz_simple import SimpleCzMatch
from masstodon.isotopes                    import isotope_calculator
from masstodon.precursor.precursor         import precursor
from masstodon.preprocessing.filters       import filter_subspectra_molecules
from masstodon.readers.from_npy            import spectrum_from_npy
from masstodon.spectrum.spectrum           import spectrum


class MasstodonBase(object):
    def __init__(self, **kwds):
        self.__dict__.update(kwds)

    def set_spectrum(self):
        self.spec = spectrum(self.mz, self.intensity)
        self.spec.bitonic_clustering()
        self.mz_digits = self.spec.bc.get_smallest_diff_digits()
        self.spec.min_mz_diff_clustering()
        self.subspectra = list(self.spec.iter_min_mz_diff_subspectra())

    def set_isotopic_calculator(self):
        self.iso_calc = isotope_calculator(digits=self.mz_digits)

    def set_molecules(self):
        self.prec = precursor(self.fasta,
                              self.charge,
                              self.name,
                              self.modifications,
                              self.fragments,
                              self.blocked_fragments,
                              self.block_prolines,
                              self.distance_charges,
                              iso_calc = self.iso_calc)
        self.mols = np.array(list(self.prec.molecules()))

    def trivial_divide_et_impera(self):
        self.good_mols, self.good_subspectra = filter_subspectra_molecules(self.subspectra,
                                                                           self.mols,
                                                                           std_cnt = self.std_cnt)

    def divide_et_impera(self):
        self.imperator = divide_ed_impera(self.good_mols,
                                          self.spec.bc,
                                          self.min_prob,
                                          self.isotopic_coverage)

    def load_imperator(self):
        self.imperator = Imperator(self.good_mols,
                                   self.spec.bc,
                                   self.min_prob,
                                   self.isotopic_coverage)
        self.imperator.load_graph(self.deconvolution_graph_path)
        self.imperator.impera()
        self.imperator.set_estimated_intensities()

    def match_estimates(self):
        self.cz_simple = SimpleCzMatch(self.good_mols, self.prec.q)
        self.cz = CzMatch(self.good_mols, self.prec.q)

    def write_csv(self, path):
        """Write results to path."""
        self.cz.write(path)
        self.cz_simple.write(path)

    def dump(self, path):
        self.spec.dump(path)
        params = {k: self.__dict__[k] for k in ('fasta',
                                                'charge',
                                                'name',
                                                'modifications',
                                                'fragments',
                                                'blocked_fragments',
                                                'block_prolines',
                                                'distance_charges',
                                                'std_cnt',
                                                'isotopic_coverage',
                                                'min_prob')}
        with open(pjoin(path, 'params.json'), 'w') as f:
            json.dump(params, f)
        self.imperator.save_graph(pjoin(path, 'deconvolution_graph.gpickle'))


def masstodon_base(mz, intensity, fasta, charge,
                   name             = "", 
                   modifications    = {},
                   fragments        = "cz",
                   blocked_fragments= ['c0'],
                   block_prolines   = True,
                   distance_charges = 5.,
                   std_cnt          = 3,
                   isotopic_coverage= .999,
                   min_prob         = .7,
                   deconvolution_graph_path = ''):
    """Run a basic session of the MassTodon.

    Parameters
    ==========
    mz : np.array
        Observed mass to charge ratios.
    intensity : np.array
        Observed intensities (corresponding to mass to charge ratios).
    fasta : str
        The FASTA sequence of the protein to study.
    charge : int
        The initial charge of the precursor filtered out in MS1.
    name : str
        The precursor's name.
    modifications : dictionary
        A dictionary of modifications.
    fragments : str
        Only 'cz' accepted for now.
        Planning other fragmentation schemes, including inner fragments.
    blocked_fragments : list
        Fragments you don't want to include, e.g. 'z5'.
    block_prolines : boolean
        Should we block prolines?
    distance_charges :
        The minimal distance between charges on the fasta sequence.
        Defaults to charges being 4 amino acids apart.
    isotopic_coverage : float
        The joint probability of the calculated isotopic distribution.
        Defaults to a decent '0.999'.
    min_prob : float
        The minimal probability an envelope has to scoop
        to be included in the deconvolution graph.
    _verbose : boolean
        Should we show the content in a verbose mode?
    """
    todon = MasstodonBase(mz = mz,
                          intensity = intensity,
                          fasta = fasta,
                          charge = charge,
                          name = name,
                          modifications = modifications,
                          fragments = fragments,
                          blocked_fragments = blocked_fragments,
                          block_prolines = block_prolines,
                          distance_charges = distance_charges,
                          std_cnt = std_cnt,
                          isotopic_coverage = isotopic_coverage,
                          min_prob = min_prob,
                          deconvolution_graph_path = deconvolution_graph_path)
    todon.set_spectrum()
    todon.set_isotopic_calculator()
    todon.set_molecules()
    todon.trivial_divide_et_impera()
    if not deconvolution_graph_path:
        todon.divide_et_impera()
    else:
        todon.load_imperator()
    todon.match_estimates()
    return todon


def masstodon_base_load(path,
                        deconvolution_graph_file = 'deconvolution_graph.gpickle'):
    mz, intensity = spectrum_from_npy(path)
    with open(pjoin(path, 'params.json'), 'r') as f:
        params = json.load(f)
    deconvolution_graph_path = pjoin(path, deconvolution_graph_file)
    params['deconvolution_graph_path'] = deconvolution_graph_path
    return masstodon_base(mz, intensity, **params)