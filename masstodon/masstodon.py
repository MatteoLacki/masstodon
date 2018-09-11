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
import numpy as np
from time import time

from masstodon.deconvolve.divide_ed_impera import divide_ed_impera, Imperator
from masstodon.estimates_matcher.cz        import CzMatch
from masstodon.estimates_matcher.cz_simple import SimpleCzMatch
from masstodon.isotopes                    import isotope_calculator
from masstodon.precursor.precursor         import precursor
from masstodon.preprocessing.filters       import filter_subspectra_molecules
from masstodon.spectrum.spectrum           import spectrum


class MasstodonBase(object):
    def set_spectrum(self, mz, intensity):
        self.spec = spectrum(mz, intensity)
        self.spec.bitonic_clustering()
        self.spec.min_mz_diff_clustering()
        self.subspectra = list(self.spec.iter_min_mz_diff_subspectra())

    def set_isotopic_calculator(self):
        mz_digits = self.spec.bc.get_smallest_diff_digits()
        self.iso_calc = isotope_calculator(digits=mz_digits)

    def set_molecules(self, fasta, charge, name, 
                      modifications, fragments,
                      blocked_fragments, block_prolines,
                      distance_charges):
        self.prec = precursor(fasta, charge, name, modifications, fragments,
                              blocked_fragments, block_prolines, distance_charges,
                              iso_calc = self.iso_calc)
        self.mols = np.array(list(self.prec.molecules()))

    def trivial_divide_et_impera(self, std_cnt=3):
        self.good_mols, self.good_subspectra = filter_subspectra_molecules(self.subspectra,
                                                                           self.mols,
                                                                           std_cnt = 3)

    def divide_et_impera(self, min_prob, isotopic_coverage):
        self.imperator = divide_ed_impera(self.good_mols,
                                          self.spec.bc,
                                          min_prob,
                                          isotopic_coverage)

    def load_imperator(self, min_prob, isotopic_coverage, deconvolution_graph_path):
        self.imperator = Imperator(self.good_mols,
                                   self.spec.bc,
                                   min_prob,
                                   isotopic_coverage)
        self.imperator.load_graph(deconvolution_graph_path)
        self.imperator.impera()
        self.imperator.set_estimated_intensities()

    def match_estimates(self):
        self.cz_simple = SimpleCzMatch(self.good_mols, self.prec.q)
        self.cz = CzMatch(self.good_mols, self.prec.q)

    def write(self, path):
        """Write results to path."""
        # self.report.write(path)
        self.cz.write(path)
        self.cz_simple.write(path)


def masstodon_base(mz, intensity, fasta, charge,
                   name             = "", 
                   modifications    = {},
                   fragments        = "cz",
                   blocked_fragments= set(['c0']),
                   block_prolines   = True,
                   distance_charges = 5.,
                   std_cnt          = 3,
                   isotopic_coverage= .999,
                   min_prob         = .7,
                   deconvolution_graph_path = '',
                   _verbose         = False):
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
    todon = MasstodonBase()
    todon.set_spectrum(mz, intensity)
    todon.set_isotopic_calculator()
    todon.set_molecules(fasta, charge, name, modifications, fragments,
                        blocked_fragments, block_prolines, distance_charges)
    todon.trivial_divide_et_impera(std_cnt)
    if not deconvolution_graph_path:
        todon.divide_et_impera(min_prob, isotopic_coverage)
    else:
        todon.load_imperator(min_prob, isotopic_coverage, deconvolution_graph_path)
    todon.match_estimates()
    return todon

