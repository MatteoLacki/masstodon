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
from collections import Counter, namedtuple
from os.path import join as pjoin

from .cz_simple       import SimpleCzMatch
from ..write.csv_tsv  import write_rows


Node = namedtuple('Node', 'type no bp q g')


class CzMatch(SimpleCzMatch):
    """Match c and z ions' intensities.

    Parameters
    ==========
    molecules : list of Molecule objects
        A list containing reaction products from one precusor.
    precursor_charge : int
        The charge of the precursor molecule.
    """
    def __init__(self,
                 molecules,
                 precursor_charge):

        self._I_ETnoD_fragments = 0
        self._I_PTR_fragments   = 0
        super().__init__(molecules, precursor_charge)

    def _get_node(self, molecule):
        """Define what should be hashed as graph node."""
        mt, po, cs = molecule._molType_position_cleavageSite()
        # TODO: this is a hack to never distinguish between HTR and ETD.
        # THIS HAS TO BE REPLACED ONCE I WILL GET RID OF TIME PRESSURE.
        if molecule.g == -1:
            g = 0
        elif molecule.g == self._Q:
            g = self._Q - 1
        else:
            g = molecule.g
        return Node(mt, po, cs, molecule.q, g)

    def _add_edge(self, C, Z):
        """Add edge between a 'c' fragment and a 'z' fragment."""
        # N_PTR = precursor.q - 1 - C.q - Z.q - C.g - Z.g
        #   N_PTR >= 0
        #       N_PTR = precursor.q - 1 - C.q - Z.q - C.g - Z.g
        #       precursor.q - 1 - C.q - Z.q - C.g - Z.g >= 0
        #       precursor.q > C.q + Z.q + C.g + Z.g
        #   N_PTR <= precursor.q - 1
        #       C.q + Z.q + C.g + Z.g >= 0  # automatically
        # N_ETnoD = C.g + Z.g
        #   N_ETnoD >= 0
        #       C.g + Z.g >= 0              # automatically
        #   N_ETnoD < precursor.q
        #       C.g + Z.g < precursor.q,
        #       but
        #       C.q + Z.q + C.g + Z.g < precursor.q
        Q = self._Q
        if C.bp == Z.bp and C.q + Z.q + C.g + Z.g < Q:
            self.graph.add_edge(C, Z, ETnoD=C.g+Z.g,
                                      PTR=Q-1-C.g-Z.g-C.q-Z.q,
                                      ETnoD_PTR=Q-1-C.q-Z.q)

    def write(self, path):
        """Write intensities and probabilities to a given path."""
        write_rows(self._iter_intensities(),
                   pjoin(path, 'pairing_intensities.csv'))
        write_rows(self._iter_probabilities(), 
                   pjoin(path, 'pairing_probabilities.csv'))

    def _get_edge_labels_4_plot(self, G):
        return {n:"${ntype}_{{{nno}}}^{{{nq}+,{ng}}}$\n{intensity}".format(
                    ntype     = n.type,
                    nno       = n.no,
                    nq        = n.q,
                    ng        = n.g,
                    intensity = d['intensity'])
                for n, d in G.nodes(data=True)}