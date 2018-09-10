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


def dict_2_string(d):
    return "".join(element + str(count) for element, count in sorted(d.items()))

def dict_2_string_readable(d):
    return " ".join(element + "_"+ str(count) for element, count in sorted(d.items()))


class Formula(LinearDict):
    pattern = get_pattern('([A-Z][a-z]?)([0-9]*)')

    @classmethod
    def recompile_pattern(cls, pattern='([A-Z][a-z]?)([0-9]*)'):
        cls.pattern = get_pattern(pattern)

    def __init__(self, formula={}):
        if isinstance(formula, str):
            formula = parse(formula, self.pattern)
        super().__init__(formula)

    def __str__(self):
        return dict_2_string_readable(self._storage)

    def __repr__(self):
        return str(self)

    def check_positivity(self):
        if any(count < 0 for element, count in self.items()):
            raise NegativeAtomCount("Attention: negative atom count after including your modifications.")

    def str_with_charges(self, q=0, g=0):
        out = self._storage.copy()
        out['H'] += q + g
        return dict_2_string(out)


def formula(formula):
    return Formula(formula)