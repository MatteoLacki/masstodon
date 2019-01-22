# -*- coding: utf-8 -*-
#   Copyright (C) 2016 Mateusz Krzysztof Łącki and Michał Startek.
#
#   This file is part of MassTodon.
#
# License: see LICENCE file.

import json
import numpy   as     np
from   os.path import join as pjoin
from   time    import time

from .data.constants import infinity
from .deconvolve.divide_ed_impera import imperator, load_imperator
from .estimates_matcher.cz import CzMatch
from .estimates_matcher.cz_simple import SimpleCzMatch
from .isotopes.calculator import isotope_calculator
from .precursor.precursor import precursor
from .preprocessing.filters import filter_subspectra_molecules
from .parse.npy import spectrum_from_npy
from .spectrum.spectrum import spectrum
from .ome.ome import ome


class Masstodon(object):
    def set_spectrum(self,
                     mz,
                     intensity,
                     min_mz_diff = 1.1,
                     orbitrap    = False,
                     threshold   = "0.0 Da"):
        self.threshold   = threshold
        self.orbitrap    = orbitrap
        self.min_mz_diff = min_mz_diff
        self.spec = spectrum(mz,
                             intensity,
                             self.min_mz_diff,
                             self.orbitrap,
                             self.threshold,
                             sort            = True,
                             drop_duplicates = True,
                             drop_zeros      = True)
        self.mz_digits = self.spec.get_smallest_diff_digits()
        self.ls        = self.spec.get_lightweight_spectrum()
        self.groups    = self.spec.get_groups()
        if min_mz_diff < infinity:
            self.subspectra = self.spec.get_min_mz_diff_subspectra()

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
        if self.min_mz_diff < infinity:
            self.good_mols, self.subspectra_within_mz_deviations = \
                self.ome.filter_by_deviations(self.subspectra,
                                              self.std_cnt)
        else:
            self.good_mols = list(self.ome.observables())

    def divide_et_impera(self, 
                         min_prob,
                         isotopic_coverage,
                         deconvolution_method = 'nnls',
                         include_zero_intensities = False):
        # store parameters
        self.min_prob                 = min_prob
        self.isotopic_coverage        = isotopic_coverage
        self.deconvolution_method     = deconvolution_method
        self.include_zero_intensities = include_zero_intensities
        # divide and rule
        self.imperator = imperator(self.good_mols,
                                   self.groups,
                                   self.ls,
                                   self.min_prob,
                                   self.isotopic_coverage,
                                   self.deconvolution_method,
                                   self.include_zero_intensities)

    def load_imperator(self, 
                       deconvolution_graph_path,
                       min_prob, 
                       isotopic_coverage,
                       deconvolution_method = 'nnls',
                       include_zero_intensities = False):
        self.deconvolution_graph_path = deconvolution_graph_path
        self.min_prob                 = min_prob
        self.isotopic_coverage        = isotopic_coverage
        self.deconvolution_method     = deconvolution_method
        self.include_zero_intensities = self.include_zero_intensities
        self.imperator = load_imperator(self.good_mols,
                                        self.groups,
                                        self.ls,
                                        self.deconvolution_graph_path,
                                        self.min_prob,
                                        self.isotopic_coverage,
                                        self.deconvolution_method,
                                        self.include_zero_intensities)

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
        if self.ome.is_one_precursor():
            self.cz.write(path)
            self.cz_simple.write(path)
        self.imperator.errors_to_json(pjoin(path, 'errors.json'))

    def dump(self, path, source_observables_graph=False, indent=4):
        """Dump the results of the fitting to locally stored files.

        Function masstodon_load can then read them back.
        """
        self.spec.dump(path) # dumping spectrum

        params = {"precursors":           self.precursors,
                  "molecules" :           self.molecules,
                  "std_cnt"   :           self.std_cnt,
                  "isotopic_coverage":    self.isotopic_coverage,
                  "min_prob"  :           self.min_prob,
                  "threshold" :           self.threshold,
                  "orbitrap"  :           self.orbitrap,
                  "min_mz_diff":          self.min_mz_diff,
                  "deconvolution_method": self.deconvolution_method,
                  "include_zero_intensities": self.include_zero_intensities}

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

    def plotly(self, *args, **kwds):
        """Save the html plotly visual output files.
        For description of the parameters, check out 
        masstodon.deconvolve.divide_ed_impera.Imperator.plotly
        """
        self.imperator.plotly(*args, **kwds)


def masstodon_batch(mz, 
                    intensity,
                    precursors          = [],
                    molecules           = [],
                    min_mz_diff         = 1.1,
                    orbitrap            = False,
                    threshold           = "0.0 Da",
                    std_cnt             = 3,
                    mz_digits           = None,
                    isotopic_coverage   = .999,
                    min_prob            = .7,
                    get_timings         = False,
                    deconvolution_method     = 'nnls',
                    include_zero_intensities = False,
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
    get_timings : boolean
        Should the output include the measured timings.
    deconvolution_graph_path : str
        A path to a valid premade deconvolution graph.
    deconvolution_method : str
        The deconvolution method to use: right now either 'nnls' or 'quantile'.
    include_zero_intensities : boolean
        Should the deconvolution include the fitting to zeros
        if some of the theoretical peaks were not found next to
        the experimental ones? Defaults to False, because it is
        not impossible that the instrument provides some intensity
        threshold for the detector.
    Returns
    =======
        An instance of the masstodon object,
        or a tuple consisting of that instance and the 
        measured timings.
    """
    t0 = time()
    m = Masstodon()
    m.set_spectrum(mz,
                   intensity,
                   min_mz_diff,
                   orbitrap,
                   threshold)
    t1 = time()
    if mz_digits is None:
        mz_digits = m.mz_digits
    m.set_isotopic_calculator()
    t2 = time()
    m.set_ome(precursors, molecules, std_cnt)
    t3 = time()
    if not deconvolution_graph_path:
        m.divide_et_impera(min_prob,
                           isotopic_coverage,
                           deconvolution_method,
                           include_zero_intensities)
    else:
        m.load_imperator(deconvolution_graph_path,
                         min_prob,
                         isotopic_coverage,
                         deconvolution_method,
                         include_zero_intensities)
    t4 = time()
    m.ome.filter_by_estimated_intensity()
    t5 = time()
    m.restrict_good_mols()
    t6 = time()
    timings = [ ("spectrum", t1-t0),
                ("isotopic_calculator", t2-t1),
                ("ome", t3-t2),
                ("imperator", t4-t3),
                ("filter_by_estimated_intensity", t5-t4),
                ("restrict_good_mols", t6-t5) ]
    if len(precursors) == 1 and len(molecules) == 0:
        m.match_estimates()
        t7 = time()
        timings.append(("match_estimates", t7-t6))
        timings.append(("total", t7-t0))
    else:
        timings.append(("total", t6-t0))
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
                     min_mz_diff       = 1.1,
                     orbitrap          = False,
                     threshold         = "0.0 Da",
                     std_cnt           = 3,
                     mz_digits         = None,
                     isotopic_coverage = .999,
                     min_prob          = .7,
                     get_timings       = False,
                     deconvolution_method     = 'nnls',
                     include_zero_intensities = False,
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
    get_timings : boolean
        Should the output include the measured timings.
    deconvolution_graph_path : str
        A path to a valid premade deconvolution graph.
    deconvolution_method : str
        The deconvolution method to use: right now either 'nnls' or 'quantile'.
    deconvolution_graph_path : str
        A path to a valid premade deconvolution graph.
    """
    precursors = [ { "fasta":               fasta,
                     "q":                   q,
                     "name":                name,
                     "modifications":       modifications,
                     "blocked_fragments":   blocked_fragments,
                     "block_prolines":      block_prolines,
                     "distance_charges":    distance_charges} ]
    return masstodon_batch(mz,
                           intensity,
                           precursors,
                           [],
                           min_mz_diff,
                           orbitrap,
                           threshold,
                           std_cnt,
                           mz_digits,
                           isotopic_coverage,
                           min_prob,
                           get_timings,
                           deconvolution_method,
                           include_zero_intensities,
                           deconvolution_graph_path)



def load_masstodon(path):
    mz, intensity = spectrum_from_npy(path)
    with open(pjoin(path, 'params.json'), 'r') as f:
        params = json.load(f)
    params['deconvolution_graph_path'] = pjoin(path,
                                               "deconvolution_graph.gpickle")
    m = masstodon_batch(mz, intensity, **params)
    return m