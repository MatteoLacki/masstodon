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
import re
from masstodon.data.elements import elements


def parser_factory(pattern='([A-Z][a-z]?)([0-9]*)'):
    """A factory to set up formula parsers.

    Parameters
    ----------
    pattern : str
        A regular expression that describes elements.

    Returns : function
        A formula parser.
    """
    formula_patern = re.compile('([A-Z][a-z]?)([0-9]*)')

    def parse(formula, pattern=formula_patern):
        """Parses chemical formula based on the class pattern definition.

        Parameters
        ----------
        formula : str
            The chemical formula string.

        Returns
        -------
        atomCnt : Counter
            A counter with elements for keys and atom counts for values.

        Examples
        --------
            >>> FP = formulaParser()
            >>> FP('C100H202')
            Counter({'C': 100, 'H': 202})
        """
        atomCnt = {}
        for elemTag, cnt in re.findall(pattern, formula):
            if not elemTag in elements:
                print("WARNING! Element tag {} ain't recognized.".format(elemTag))
            if cnt == '':
                cnt = 1
            else:
                cnt = int(cnt)
            atomCnt[elemTag] = cnt
        return atomCnt

    return parse

parser = parser_factory()