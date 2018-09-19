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

from masstodon.deconvolve.divide_ed_impera import imperator, load_imperator
from masstodon.estimates_matcher.cz        import CzMatch
from masstodon.estimates_matcher.cz_simple import SimpleCzMatch
from masstodon.isotopes                    import isotope_calculator
from masstodon.precursor.precursor         import precursor
from masstodon.preprocessing.filters       import filter_subspectra_molecules
from masstodon.readers.from_npy            import spectrum_from_npy
from masstodon.spectrum.spectrum           import spectrum
from masstodon.ome.ome                     import ome

class Masstodon(object):
    def set_spectrum(self, mz, intensity):
        self.spec = spectrum(mz, intensity)
        self.spec.bitonic_clustering()
        self.mz_digits = self.spec.bc.get_smallest_diff_digits()
        self.spec.min_mz_diff_clustering()
        self.subspectra = list(self.spec.iter_min_mz_diff_subspectra())

    def set_isotopic_calculator(self):
        self.iso_calc = isotope_calculator(digits=self.mz_digits)

    def set_ome(self, 
                precursors=[],
                molecules=[],
                std_cnt=3):
        self.precursors = precursors
        self.molecules  = molecules
        self.std_cnt    = std_cnt
        self.ome = ome(self.iso_calc,
                       self.precursors,
                       self.molecules)
        self.good_mols, self.good_subspectra = \
            self.ome.filter_by_deviations(self.subspectra,
                                          self.std_cnt)

    def divide_et_impera(self, 
                         min_prob,
                         isotopic_coverage):
        self.min_prob          = min_prob
        self.isotopic_coverage = isotopic_coverage
        self.imperator = imperator(self.good_mols,
                                   self.spec.bc,
                                   self.min_prob,
                                   self.isotopic_coverage)

    def load_imperator(self, 
                       deconvolution_graph_path,
                       min_prob, 
                       isotopic_coverage):
        self.deconvolution_graph_path = deconvolution_graph_path
        self.min_prob                 = min_prob
        self.isotopic_coverage        = isotopic_coverage
        self.imperator = load_imperator(self.good_mols,
                                        self.spec.bc,
                                        self.deconvolution_graph_path,
                                        self.min_prob,
                                        self.isotopic_coverage)

    def match_estimates(self):
        """Match the fragment intensities to get the idea about 
        the minimal number of precursor molecules."""
        sources = self.ome.sources()
        prec = next(sources)
        try:
            x = next(sources)
            raise AttributeError("You supplied too many precursors for the c/z analysis.")
        except StopIteration:
            pass
        for mol in self.ome.observables():
            mol.name = self.ome.G[mol][prec]['name']
            mol.prec_fasta_len = len(prec.fasta)
        self.cz_simple = SimpleCzMatch(self.good_mols, prec.q)
        self.cz = CzMatch(self.good_mols, prec.q)

    def write(self, path):
        """Write results to path."""
        self.ome.write(pjoin(path, 'estimates.csv'))
        self.cz.write(path)
        self.cz_simple.write(path)
        self.imperator.errors_to_json(pjoin(path, 'errors.json'))

    def dump(self, path, source_observables_graph=False, indent=None):
        """Dump the results of the fitting to locally stored files.

        Function masstodon_load can then read them back.
        """
        self.spec.dump(path)
        params = {"precursors": self.precursors,
                  "molecules" : self.molecules,
                  "std_cnt"   : self.std_cnt,
                  "isotopic_coverage" : self.isotopic_coverage,
                  "min_prob"  : self.min_prob}
        with open(pjoin(path, 'params.json'), 'w') as f:
            json.dump(params, f, indent=indent)
        self.imperator.save_graph(pjoin(path, 'deconvolution_graph.gpickle'))
        self.imperator.errors_to_json(pjoin(path, 'errors.json'),
                                      indent=indent)
        if source_observables_graph:
            self.ome.dump(pjoin(path,
                                'sources_observables_graph.gpickle'))
        self.ome.dump_stats(pjoin(path, 'deconvolution_graph_stats.json'),
                            indent=indent)

    def restrict_good_mols(self):
        self.good_mols = [m for m in self.good_mols if m.intensity > 0.0]

    def plotly(self, path, show=True):
        """Save the html plotly visual output files.

        Paramaters
        ==========
        path : str
            Output folder.
        """
        self.imperator.plotly_solutions(pjoin(path,
                                             "spectrum.html"),
                                        show)


def masstodon_batch(mz, 
                    intensity, 
                    precursors          = [],
                    molecules           = [],
                    std_cnt             = 3,
                    mz_digits           = None,
                    isotopic_coverage   = .999,
                    min_prob            = .7,
                    get_timings         = False,
                    deconvolution_graph_path = ''):
    """Run a session of masstodon with multiple sources.

    Parameters
    ==========
    mz : np.array
        Observed mass to charge ratios.
    intensity : np.array
        Observed intensities (corresponding to mass to charge ratios).
    precursors : list of dicts
        Arguments to the precursor function.
    molecules : list of dicts
        Arguments to the molecule function.
    all_molecules_kwds : dict
        Arguments
    std_cnt : float
        Number of standard deviations around average theoretical m/z of a fragment to be considered in the trivial filtering.
    mz_digits : int
        Number of significant digits while rounding m/z values.
    isotopic_coverage : float
        The joint probability of the calculated isotopic distribution.
        Defaults to a decent '0.999'.
    min_prob : float
        The minimal probability an envelope has to scoop
        to be included in the deconvolution graph.
    deconvolution_graph_path : str
        A path to a valid premade deconvolution graph.
    """
    t0 = time()
    m = Masstodon()
    m.set_spectrum(mz, intensity)
    t1 = time()
    if mz_digits is None:
        mz_digits = m.mz_digits
    m.set_isotopic_calculator()
    t2 = time()
    m.set_ome(precursors, molecules, std_cnt)
    t3 = time()
    if not deconvolution_graph_path:
        m.divide_et_impera(min_prob, isotopic_coverage)
    else:
        m.load_imperator(deconvolution_graph_path,
                         min_prob,
                         isotopic_coverage)
    t4 = time()
    m.ome.filter_by_estimated_intensity()
    t5 = time()
    m.restrict_good_mols()
    t6 = time()
    if m.ome.is_one_precursor():
        m.match_estimates()
        t7 = time()
    timings = (
        ("spectrum", t1-t0),
        ("isotopic_calculator", t2-t1),
        ("ome", t3-t2),
        ("imperator", t4-t3),
        ("filter_by_estimated_intensity", t5-t4),
        ("restrict_good_mols", t6-t5),
        ("match_estimates", t7-t6),
        ("total", t7-t0)
    )
    if get_timings:
        return m, timings
    else:
        return m


def masstodon_single(mz, intensity, fasta, q,
                     name              = "", 
                     modifications     = {},
                     fragments         = "cz",
                     blocked_fragments = ['c0'],
                     block_prolines    = True,
                     distance_charges  = 5.,
                     std_cnt           = 3,
                     mz_digits         = None,
                     isotopic_coverage = .999,
                     min_prob          = .7,
                     get_timings       = False,
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
    q : int
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
    std_cnt : float
        Number of standard deviations around average theoretical m/z of a fragment to be considered in the trivial filtering.
    mz_digits : int
        Number of significant digits while rounding m/z values.
    isotopic_coverage : float
        The joint probability of the calculated isotopic distribution.
        Defaults to a decent '0.999'.
    min_prob : float
        The minimal probability an envelope has to scoop
        to be included in the deconvolution graph.
    deconvolution_graph_path : str
        A path to a valid premade deconvolution graph.
    """
    precursors = [ { "fasta":               fasta,
                     "q":                   q,
                     "name":                name,
                     "modifications":       modifications,
                     "blocked_fragments":   blocked_fragments,
                     "block_prolines":      block_prolines,
                     "distance_charges":    distance_charges } ]
    return masstodon_batch(mz, intensity, precursors,
                           std_cnt                  = std_cnt,
                           mz_digits                = mz_digits,
                           isotopic_coverage        = isotopic_coverage,
                           min_prob                 = min_prob,
                           get_timings              = get_timings,
                           deconvolution_graph_path = deconvolution_graph_path)



def masstodon_load(path):
    mz, intensity = spectrum_from_npy(path)
    with open(pjoin(path, 'params.json'), 'r') as f:
        params = json.load(f)
    params['deconvolution_graph_path'] = pjoin(
        path, "deconvolution_graph.gpickle")
    m = masstodon_batch(mz, intensity, **params)
    return m