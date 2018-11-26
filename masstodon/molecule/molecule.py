# -*- coding: utf-8 -*-
#
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
from masstodon.data.constants  import infinity
from masstodon.plot.spectrum   import plot_spectrum
from masstodon.formula.formula import as_formula


class Molecule(object):
    def __init__(self, formula, iso_calc, q=0, g=0):
        self.formula   = as_formula(formula)
        self.q         = int(q)
        self.g         = int(g)
        self.intensity = 0.0
        self.iso_calc  = iso_calc

    # TODO generalize to abxy
    def _molType_position_cleavageSite(self):
        """Supply information necessary for the matching of the estimated intensities.

        Parameters
        ==========
        mol_name : str
            The name of the molecule as induced by the precursor.
        Returns
        =======
        tuple : type of molecule (p-recursor, c-fragment, z-fragment),
                position in the fasta file,
                number of the c or z fragment or None for precursor.
        """
        mt = self.name[0]
        if mt == 'p':
            return None
        else:
            po = int(self.name[1:])
            cs = None if mt == 'p' else \
                   po if mt == 'c' else self.prec_fasta_len - po
            return mt, po, cs

    @property
    def monoisotopic_mz(self):
        return self.iso_calc.monoisotopic_mz(self.formula,
                                             self.q,
                                             self.g)

    @property
    def mean_mz(self):
        return self.iso_calc.mean_mz(self.formula,
                                     self.q,
                                     self.g)

    @property
    def sd_mz(self):
        return self.iso_calc.sd_mz(self.formula,
                                   self.q,
                                   self.g)

    def interval(self, std_cnt = 4):
        mean_mz = self.mean_mz
        sd_mz   = self.sd_mz
        s = mean_mz - std_cnt * sd_mz
        e = mean_mz + std_cnt * sd_mz
        return s, e

    def isotopologues(self,
                      prob    = .999,
                     _memoize = True):
        return self.iso_calc(self.formula,
                             prob,
                             self.q,
                             self.g,
                            _memoize=_memoize)

    def __repr__(self):
        return "({f} q={q} g={g} I={I_int})".format(
            I_int = int(self.intensity),
            f = self.formula.str_with_charges(self.q, self.g),
            **self.__dict__)

    def __hash__(self):
        """The least you need to know to trace a molecule.

        The molecule is uniquely defined by its total atom count and charge.
        As best summarized by Metallica: nothing else matters.
        """
        return hash((self.formula.str_with_charges(self.q, self.g),
                     self.q))

    def __eq__(self, other):
        A = self.formula.str_with_charges(self.q, self.g) == \
            other.formula.str_with_charges(other.q, other.g)
        B = self.q == other.q
        return A and B 

    def plot(self,
             plt_style = 'dark_background',
             show      = True):
        """Plot the molecules isotopic distribution."""
        env = self.isotopologues()
        env.plot(plt_style = 'dark_background',
                 show      = show)


def molecule(formula, iso_calc, q=0, g=0):
    mol = Molecule(formula, iso_calc, q, g)
    return mol

