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
from cvxopt import  matrix, spmatrix, sparse, spdiag, solvers
import networkx as nx
from networkx import connected_component_subgraphs as connected_components

from masstodon.data.constants import eps, infinity
from masstodon.write.csv_tsv  import write_rows


def diag(val, dim):
    """Make a sparse identity matrix multiplied by a scalar val."""
    return spdiag([spmatrix(val,[0],[0]) for i in range(dim)])


def incidence_matrix(graph, row_cnt, col_cnt):
    """Make a sparse incidence matrix of the graph G."""
    L = spmatrix([], [], [], size=(row_cnt, col_cnt) )
    NodesNo = dict([ (N,i) for i,N in enumerate(graph)])
    for j, (N0, N1) in enumerate(graph.edges()):
        L[NodesNo[N0],j] = 1
        L[NodesNo[N1],j] = 1
    return L


Node = namedtuple('Node', 'type no bp q')

class SimpleCzMatch(object):
    """Match c and z ions' intensities neglecting the quenched charge.

    Parameters
    ==========
    molecules : list of Molecule objects
        A list containing reaction products from one precusor.
    precursor_charge : int
        The charge of the precursor molecule.
    show_progress : boolean
        Show progress of the CVXOPT calculations.
    maxiters : int
        Maximum number of iterations for the CVXOPT algorithm.

    """
    def __init__(self,
                 molecules,
                 precursor_charge,
                 show_progress=False,
                 maxiters=1000,
                 **kwds):
        solvers.options['show_progress'] = bool(show_progress)
        solvers.options['maxiters'] = int(maxiters)

        self._molecules = molecules
        self._Q = int(precursor_charge)
        # _I_ = Intensity
        self._I_ETDorHTR_bond = Counter()
        self._I_ETDorHTR = 0.0
        self._I_ETnoD_precursor = 0.0
        self._I_PTR_precursor = 0.0
        self._I_ETnoD_PTR_precursor = Counter()  # len(ETnoD), len(PTR) -> Intensity
        self._I_ETnoD_PTR_fragments = 0.0
        self._I_ETnoD_PTR_bond = Counter()
        self._I_lavish = 0.0
        self._make_graph()
        # self._match()
        # self.intensities = self._get_intensities()
        # self.probabilities = self._get_probabilities()
        # self.branching_ratio = self._I_PTR_precursor / self._I_ETnoD_precursor if self._I_ETnoD_precursor else None

    def _get_node(self, molecule):
        """Define what should be hashed as graph node."""
        mt, po, cs = molecule._molType_position_cleavageSite()
        return Node(mt, po, cs, molecule.q)

    def _add_edge(self, C, Z):
        """Add edge between a 'c' fragment and a 'z' fragment."""
        Q = self._Q
        if C.bp == Z.bp and C.q + Z.q < Q:
            self.graph.add_edge(C, Z, ETnoD_PTR = Q - 1 - C.q - Z.q)

    def _add_self_loop(self, N):
        """Add edge between a 'c' fragment and a 'z' fragment."""
        self.graph.add_edge(N, N, ETnoD_PTR=self._Q - 1 - N.q)

    def _make_graph(self):
        """Prepare the matching graph."""
        Q = self._Q
        self.graph = nx.Graph()
        for mol in self._molecules:
            estimate = mol.intensity
            if estimate > 0:
                if mol.name == 'precursor':
                    g = mol.g
                    q = mol.q
                    self._I_ETnoD_precursor += g * estimate
                    self._I_PTR_precursor += (Q - q - g) * estimate
                    self._I_ETnoD_PTR_precursor[g, Q - q - g] = estimate
                else:
                    frag = self._get_node(mol)
                    if not frag in self.graph:
                        self.graph.add_node(frag, intensity=0)
                        # the next line must be one indentation lower!!!
                        # don't you forget about it!!!
                    self.graph.node[frag]['intensity'] += estimate
        for C in self.graph:
            if C.type == 'c':
                for Z in self.graph:
                    if Z.type == 'z':
                        self._add_edge(C, Z)
        for N in self.graph:
            self._add_self_loop(N)

    def _optimize(self, G):
        """Match the intensities in a cluster.

        Decorates the self.graph with flow values.
        """
        Q = self._Q
        # lavish: all fragments lose cofragments
        lavish = sum((Q - 1 - N.q) * I for N, I in G.nodes.data('intensity'))
        self._I_lavish += lavish
        if len(G) > 1:
            intensities = matrix([float(I) for N, I in G.nodes.data('intensity')])
            costs = matrix([float(ETnoD_PTR) for N,M, ETnoD_PTR
                            in G.edges.data('ETnoD_PTR')])
            edges_cnt    = G.size()  # number of c-z pairings
            equalities   = incidence_matrix(G, len(intensities), edges_cnt)
            inequalities = diag(-1.0, edges_cnt)
            upper_bounds = matrix([0.0] * edges_cnt)
            primalstart = {}
            primalstart['x'] = matrix([0.0] * edges_cnt)
            primalstart['s'] = matrix([eps] * len(upper_bounds))
            solution = solvers.conelp(c = costs,
                                      G = inequalities,
                                      h = upper_bounds,
                                      A = equalities,
                                      b = intensities,
                                      primalstart = primalstart)
            self._I_ETnoD_PTR_fragments += solution['primal objective']
            for i, (N, M) in enumerate(G.edges()):
                self.graph[N][M]['flow'] = solution['x'][i]
        else:
            self._I_ETnoD_PTR_fragments += lavish
            N, N_intensity = list(G.nodes.data('intensity'))[0]
            self.graph[N][N]['flow'] = N_intensity

    def _get_intensities(self):
        """Estimate intensities."""
        # _I_ = Intensity
        assert self._I_ETnoD_PTR_fragments >= 0,\
            "The total intensity of ETnoD and PTR on fragments should be non-negative.\
            But it equals {}".format(self._I_ETnoD_PTR_fragments)

        for N, M, data in self.graph.edges.data():
            self._I_ETDorHTR_bond[M.bp] += data['flow'] # intensity
            if data['ETnoD_PTR']:
                ETnoD_or_PTR_fragment = data['ETnoD_PTR'] * data['flow']
                self._I_ETnoD_PTR_bond[M.bp] += ETnoD_or_PTR_fragment

        self._I_ETDorHTR = sum(v for k, v in self._I_ETDorHTR_bond.items())
        self._I_total_ETnoDorPTR = self._I_ETnoD_precursor + \
                                   self._I_PTR_precursor + \
                                   self._I_ETnoD_PTR_fragments
        self._I_reactions = self._I_total_ETnoDorPTR + \
                            self._I_ETDorHTR
        self._I_unreacted_precursor = self._I_ETnoD_PTR_precursor[0,0]

        return {k[3:]: v for k, v in self.__dict__.items() if k[0:3] == '_I_'}

    def _iter_intensities(self):
        """Generate rows for a csv/tsv file with estimated intensities."""
        yield ('unreacted', 'total:', self._I_unreacted_precursor)
        yield ('reacted', 'total:', self._I_reactions)
        yield ('ETD or HTR', 'total:', self._I_ETDorHTR)
        bonds = list(self._I_ETDorHTR_bond.items())
        bonds.sort()
        for no, v in bonds:
            yield ('', 'bond %d' % no, v)
        yield ('ETnoD or PTR on fragments', 'total:', self._I_ETnoD_PTR_fragments)
        bonds = list(self._I_ETnoD_PTR_bond.items())
        bonds.sort()
        for no, v in bonds:
            yield ('', 'bond %d' % no, v)
        yield ('ETnoD on precursors', 'total:', self._I_ETnoD_precursor)
        yield ('PTR on precursors', 'total:', self._I_PTR_precursor)
        yield ('ETnoD or PTR', 'total:', self._I_total_ETnoDorPTR)

    def _get_probabilities(self):
        """Estimate probabilities."""
        # _P_ = Probability
        if self._I_ETDorHTR > 0:
            self._P_fragmentation_bond = {k: v/self._I_ETDorHTR for k, v in \
                                          self._I_ETDorHTR_bond.items()}

        if self._I_reactions + self._I_unreacted_precursor > 0:
            self._P_reaction = self._I_reactions / (self._I_reactions + self._I_unreacted_precursor)

        if self._I_reactions > 0:
            self._P_fragmentation = self._I_ETDorHTR / self._I_reactions

        ETnoD_prec = float(self._I_ETnoD_precursor)
        PTR_prec = float(self._I_PTR_precursor)
        if ETnoD_prec + PTR_prec > 0:
            self._P_ETnoD_precursor = ETnoD_prec / (ETnoD_prec + PTR_prec)
            self._P_PTR_precursor = 1.0 - self._P_ETnoD_precursor

        if self._I_total_ETnoDorPTR > 0:
            self._P_ETnoD_PTR = self._I_total_ETnoDorPTR / self._I_reactions

        return {k[3:]: v for k, v in self.__dict__.items() if k[0:3] == '_P_'}

    def _iter_probabilities(self):
        """Generate rows for a csv/tsv file with estimated probabilities."""

        if self._I_ETnoD_precursor + self._I_PTR_precursor > 0:
            yield ('probability', 'ETnoD (precusor)', "{:10.3f}%".format(100 * self._P_ETnoD_precursor))
            yield ('probability', 'PTR (precusor)', "{:10.3f}%".format(100 * self._P_PTR_precursor))

        if self._I_total_ETnoDorPTR > 0:
            yield ('probability', 'ETnoD of PTR', "{:10.3f}%".format(100 * self._P_ETnoD_PTR))

        if self._I_reactions > 0:
            yield ('probability', 'fragmentation', "{:10.3f}%".format(100 * self._P_fragmentation))
            bonds = list(self._P_fragmentation_bond.items())
            bonds.sort()
            for no, v in bonds:
                yield ('', 'bond %d' % no, "{:10.3f}%".format(100 * v) )

    def write(self, path):
        """Write intensities and probabilities to a given path."""
        write_rows(self._iter_intensities(), path + 'simple_pairing_intensities.csv')
        write_rows(self._iter_probabilities(), path + 'simple_pairing_probabilities.csv')

    def _match(self):
        """Pair molecules minimizing the number of reactions and calculate the resulting probabilities."""
        for G in connected_components(self.graph):
            self._optimize(G)
