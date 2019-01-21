# -*- coding: utf-8 -*-
#   Copyright (C) 2016 Mateusz Krzysztof Łącki and Michał Startek.
#
#   This file is part of MassTodon.
#
# License: see LICENCE file.

import re
from collections import defaultdict

# check it out on https://www.debuggex.com/
chem_diff = re.compile(r"([A-Z][a-z]?)(-?\d*)")
mod_pattern = re.compile(r"(\d+)(C_carbo|N|C_alpha)?")


def parse(mod):
    """Parse a modification to the amino acids.

    Parameters
    ==========
    mod : str
        The modification consists of its position and chemical diff separated by '=' sign, e.g. "11=H-1N1O-1". The position can be further limited to a group in the back-chain of the amino acid (by default, the modification changes only the side-chain, i.e. changes the group of the C_alpha). Other possible atoms are denoted by 'N' and 'C_carbo'.

    Returns
    =======
    Triplet (position, atom, chemical diff), i.e. (11, 'C_carbo', {'H': 1, 'N': 1, 'O': -1}).
    """
    pos, m = mod.split("=")
    pos = re.match(mod_pattern, pos)

    aa_no = int(pos[1])
    aa_atom = pos[2] if pos[2] else "C_alpha"

    diff = {a: 1 if c is '' else int(c) for a, c in re.findall(chem_diff, m)}

    return aa_no, aa_atom, diff


def test_parse():
    assert parse("11=C-1H20Ag10") == (11, 'C_alpha', {'C':-1, 'H':20, 'Ag':10})
    assert parse("10N=CAr10") == (10, 'N', {'C':1, 'Ar':10})
    assert parse("10C_carbo=P10") == (10, 'C_carbo', {'P':10})
    assert parse("10C_alpha=H-20P10") == (10, 'C_alpha', {'H':-20, 'P':10})


def parse_mods(mods):
    out = defaultdict(lambda:defaultdict(dict))
    for mod in mods:
        aa_no, aa_atom, diff = parse(mod)
        out[aa_no][aa_atom] = diff
    out = dict(out)
    for i in out:
        out[i] = dict(out[i])
    return out
