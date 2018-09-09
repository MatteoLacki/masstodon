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
from collections    import Counter
from math           import ceil, log10
from time           import time
import csv
import os

# from masstodon.data.constants             import eps
# from masstodon.ceconvolution.Deconvolve   import build_deconvolution_graph
# from masstodon.Deconvolution.Deconvolve   import deconvolve
# from masstodon.MatchMaker.CzMatch         import CzMatch
# from masstodon.MatchMaker.SimpleCzMatch   import SimpleCzMatch
# from masstodon.Parsers.Paths              import parse_path
# from masstodon.Plot.bokeh_spectrum        import bokeh_spectrum
# from masstodon.Precursor.Precursor        import Precursor
# from masstodon.Spectra.Spectrum           import Spectrum
# from masstodon.Reporter.Reporter          import Reporter

# #TODO: simplify this call: there are way too many explicit arguments
# class MassTodon(object):
#     def __init__(self,
#                  spectrum,
#                  fasta,
#                  charge,
#                  mz_tol,
#                  mz_digits=-10.,
#                  name="",
#                  modifications={},
#                  fragments="cz",
#                  blocked_fragments=set(['c0']),
#                  block_prolines=True,
#                  distance_charges=5.,
#                  min_intensity=eps,
#                  percent_top_peaks=.999,
#                  deconvolution_method='Matteo',
#                  joint_probability=.999,
#                  min_prob_per_molecule=.7,
#                  _max_buffer_len=0.5,
#                  _L1_flow=0.01,
#                  _L2_flow=0.01,
#                  _L1_intensity=0.01,
#                  _L2_intensity=0.01,
#                  _max_times=10.,
#                  _show_progress=False,
#                  _maxiters=1000.,
#                  _verbose=False,
#                  sigma2=.1,
#                  ni2=.1):
#         """Run a full session of the MassTodon.

#         Parameters
#         ==========
#         spectrum : st, tuple of numpy.arrays, or an instance of Spectrum
#             The string path can end with:
#                 *.txt, for spectra saved in tab separated format.
#                 *.mzXml, for spectra saved in mzXml format.
#                 *.mzml, for spectra saved mzml format.
#             The tuple consists of two numpy arrays:
#                 one with m over z ratios,
#                 the other with intensities.
#         fasta : str
#             The FASTA sequence of the protein to study.
#         charge : int
#             The initial charge of the precursor filtered out in MS1.
#         mz_tol : float
#             The tolerance in the m/z axis.
#             Ultimately turned into a function mz_tol(mz),
#             reporting the tolerance interval - a tuple '(mz_L, mz_R)',
#             where mz_L = mz - mz_tol and mz_R = mz + mz_tol.
#             Accepts user-provided functions too: in that case, the header should be
#             'def mz_tol(mz):' and it should return a tuple '(mz_L, mz_R)'.
#             If providing a function, it is necessary to provide
#             mz_digits.
#         mz_digits : int
#             The number of significant digits used for m/z representation.
#             If not explicitly provided, set to 'int(ceil(-log10(mz_tol)))'.
#             Defaults to nonsensical -10.
#         name : str
#             The precursor's name.
#         modifications : dictionary
#             A dictionary of modifications.
#         fragments : str
#             Only 'cz' accepted for now.
#             Planning other fragmentation schemes, including inner fragments.
#         blocked_fragments : list
#             Fragments you don't want to include, e.g. 'z5'.
#         block_prolines : boolean
#             Should we block prolines?
#         distance_charges :
#             The minimal distance between charges on the fasta sequence.
#             Defaults to charges being 4 amino acids apart.
#         min_intensity : float
#             Experimental peaks with lower height will be trimmed.
#         percent_top_peaks : float
#             Percentage of the heighest peaks in the spectrum to be included.
#         deconvolution_method : str
#             Input 'Matteo' for MassTodon paper deconvolution.
#             Input 'Ciacho_Wanda' for gaussian kernel deconvolution.
#         joint_probability : float
#             The joint probability of the calculated isotopic distribution.
#             Defaults to a decent '0.999'.
#         min_prob_per_molecule : float
#             The minimal probability an envelope has to scoop
#             to be included in the deconvolution graph.
#         _max_buffer_len : float
#             The maximal length of the visual buffer between peaks, i.e.
#             the big rectangle width.
#         _L1_flow : float
#             L1 penalty for high flows of intensities.
#         _L2_flow : float
#             L2 penalty (a.k.a. ridge regression like) for high flows of intensities.
#         _L1_intensity : float
#             L1 penalty for high intensity estimates.
#         _L2_intensities : float
#             L2 penalty (a.k.a. ridge regression like) for high intensities.
#         _max_times : int
#             The maximal number of times to run CVXOPT.
#         _show_progress : boolean
#             Show progress of the CVXOPT calculations.
#         _maxiters : int
#             Maximum number of iterations for the CVXOPT algorithm.
#         _verbose : boolean
#             Should we show the content in a verbose mode?
#         sigma2 : float
#             Variance of the experimental peak's m/z ratio.
#         ni2 : float
#             Variance of the theoretic isotopologue's m/z ratio.
#         """
#         self.times = {}
#         t0 = time()
#         if isinstance(mz_tol, str):
#             mz_tol = float(mz_tol)
#         if not isinstance(mz_tol, float):
#             assert isinstance(mz_digits, int) and mz_digits != -10
#         self.mz_digits = int(ceil(-log10(mz_tol))) if mz_digits in (-10,'-10') else mz_digits
#         self.precursor = Precursor(name=name,
#                                    fasta=fasta,
#                                    charge=charge,
#                                    modifications=modifications,
#                                    fragments=fragments,
#                                    blocked_fragments=blocked_fragments,
#                                    block_prolines=block_prolines,
#                                    distance_charges=distance_charges)
#         self.molecules = list(self.precursor.molecules())
#         t1 = time()
#         self.times["creating_molecules"] = t1 - t0
#         self.spectrum = Spectrum(spectrum=spectrum,
#                                  mz_digits=self.mz_digits,
#                                  min_intensity=min_intensity,
#                                  percent_top_peaks=percent_top_peaks)
#         t2 = time()
#         self.times["spectrum_preprocessing"] = t2 - t1
#         self.deconvolution_method = deconvolution_method
#         self.deconvolution_graph = \
#             build_deconvolution_graph(molecules=self.molecules, 
#                                       spectrum=self.spectrum,
#                                       mz_tol=mz_tol,
#                                       min_prob_per_molecule=min_prob_per_molecule,
#                                      _verbose=_verbose,
#                                       method=deconvolution_method,
#                                       # IsoSpec arguments:
#                                       joint_probability=joint_probability, 
#                                       mz_digits=self.mz_digits)
#         t3 = time()
#         self.times["deconvolution_graph"] = t3 - t2
#         self.solutions = \
#             deconvolve(deconvolution_graph=self.deconvolution_graph,
#                        method=self.deconvolution_method,
#                       _verbose=_verbose,
#                        # deconvutor arguments
#                        L1_flow=_L1_flow,
#                        L2_flow=_L2_flow,
#                        L1_intensity=_L1_intensity,
#                        L2_intensity=_L2_intensity,
#                        max_times=_max_times,
#                        show_progress=_show_progress,
#                        maxiters=_maxiters)
#         #TODO: leaving as generator causes problems: no 'len' to call later on.
#         self.solutions = list(self.solutions)
#         t4 = time()
#         self.times["deconvolution"] = t4 - t3
#         self.report = Reporter(masstodon=self,
#                                max_buffer_len=_max_buffer_len)
#         t5 = time()
#         self.times["report"] = t5 - t4
#         #TODO: change the code below so that it could handle
#         #      reaction products from different precursors
#         #TODO: make it work with a,b,x,y fragment types
#         self.simple_cz_match = \
#             SimpleCzMatch(molecules=self.molecules,
#                           precursor_charge=self.precursor.q)
#         t6 = time()
#         self.times["simple_cz_match"] = t6 - t5
#         self.cz_match = \
#             CzMatch(molecules=self.molecules,
#                     precursor_charge=self.precursor.q)
#         t7 = time()
#         self.times["cz_match"] = t7 - t6


#     def write(self, path):
#         """Write results to path."""
#         self.report.write(path)
#         self.simple_cz_match.write(path)
#         self.cz_match.write(path)
#         with open(os.path.join(path, 'dict.csv'), 'w') as h:
#             writer = csv.writer(h)
#             for key, value in self.times.items():
#                 writer.writerow([key, value])

