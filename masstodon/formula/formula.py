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
from masstodon.formula.parse import get_pattern, parse
from masstodon.formula.linear_dict import LinearDict


class NegativeAtomCount(Exception):
    pass


def dict2string(d):
    return "".join(element + str(count) for element, count in sorted(d.items()))


def dict2string_readable(d):
    return " ".join(element + "_"+ str(count) for element, count in sorted(d.items()))


def dict2tex(d):
    f = "".join(f"{element}_{{{count}}}" if count > 1 else element
                for element, count in sorted(d.items()))
    return f"${f}$"


class Formula(LinearDict):
    pattern = get_pattern('([A-Z][a-z]?)([0-9]*)')

    @classmethod
    def recompile_pattern(cls, pattern='([A-Z][a-z]?)([0-9]*)'):
        cls.pattern = get_pattern(pattern)

    def __init__(self, formula={}):
        """Initialize the formula.

        Parameters
        ==========
        formula : str or dict
            Either a string, e.g. 'C100H202', or a dictionary
            specifying the chemical formula, e.g. "{'C':100, 'H':202}".
        """
        if isinstance(formula, str):
            formula = parse(formula, self.pattern)
        super().__init__(formula)

    def __str__(self):
        return dict2string(self._storage)

    def __repr__(self):
        return str(self)

    def check_positivity(self):
        if any(count < 0 for element, count in self.items()):
            raise NegativeAtomCount("Attention: negative atom count after including your modifications.")

    def str_with_charges(self, q=0, g=0):
        """String that includes charges and quenched charges.

        This is mainly used to uniquely hash a molecule.
        """
        out = self._storage.copy()
        out['H'] += q + g
        return dict2string(out)

    def tex_with_charges(self, q=0, g=0):
        out = self._storage.copy()
        out['H'] += q + g
        return dict2tex(out)

    def tex(self):
        return dict2tex(out)


def as_formula(f):
    if isinstance(f, str):
        return Formula(f)
    else:
        return f