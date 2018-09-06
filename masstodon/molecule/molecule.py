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
from masstodon.data.constants import infinity
from masstodon.plot.spectrum  import plot_spectrum
from masstodon.isotopes       import iso_calc


class Molecule(object):
    def __init__(self, name,
                       source,
                       formula,
                       iso_calc = iso_calc,
                       q        = 0,
                       g        = 0):
        self.name      = name
        self.source    = source
        self.formula   = formula
        self.q         = int(q)
        self.g         = int(g)
        self.intensity = 0.0
        self.iso_calc  = iso_calc

    #TODO generalize to abxy
    def _molType_position_cleavageSite(self):
        mt = self.name[0]
        if mt == 'p':
            return None
        else:
            po = int(self.name[1:])
            fasta_len = len(self.source.fasta)
            cs = None if mt == 'p' else \
                   po if mt == 'c' else fasta_len - po
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
        return "({name} {source.name} q={q} g={g} I={I_int})".format(
            I_int=int(self.intensity),
            **self.__dict__)

    def __hash__(self):
        return hash((self.name,
                     self.source,
                     str(self.formula),
                     self.q,
                     self.g))

    def plot(self,
             plt_style = 'dark_background',
             show      = True):
        """Plot the molecules isotopic distribution."""
        env = self.isotopologues()
        env.plot(plt_style = 'dark_background',
                 show      = show)


def molecule(name,
             source,
             formula,
             iso_calc = iso_calc,
             q        = 0,
             g        = 0):
    mol = Molecule(name, source, formula, iso_calc, q, g)
    return mol

