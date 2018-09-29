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
from masstodon.data.amino_acids     import get_amino_acids
from masstodon.formula.formula      import Formula
from masstodon.molecule.molecule    import Molecule


def flatten_modification(mod):
    return ''.join(''.join((str(n), t, str(v)))
                   for (n, t), v in mod.items())

class Precursor(Molecule):
    """Make precursor.

    Parameters
    ==========
    fasta : str
        The fasta of the studied molecular species.
    name : str
        The name of the precursor molecule.
    charge : int
        The charge of the precursor ion.
    modifications : dictionary
        A dictionary of modifications of amino acids,
        with key equal to amino acid number in fasta sequence
        (from N to C termini), and value equal to the dictionary 
        with group modifications.
        The modifications have keys equalt to C_carbo, C_alpha, or N, 
        and values equal to atom count in form of a linearCounter.
    fragments: str
        For now 'cz' only, but we are working on it.
    blocked_fragments : list
        Fragments you don't want to include, e.g. 'z5'.
    block_prolines : boolean
        Should we block prolines?
    distance_charges :
        The minimal distance between charges on the fasta sequence.
        Defaults to charges being 4 amino acids apart.
    kwds :
        Settings for other methods.
    """
    amino_acids = get_amino_acids()  # residues only!

    def __init__(self,
                 fasta,
                 q,
                 iso_calc,
                 name              = "",
                 modifications     = {},
                 fragments         = "cz",
                 blocked_fragments = set(['c0']),
                 block_prolines    = True,
                 distance_charges  = 5,
                 **kwds):
        self.name = name
        self.real = True
        self.fasta = fasta
        self.q     = int(q)
        self.g     = 0 # To get the plotting from superclass.
        self.iso_calc = iso_calc
        self.groups   = ('N', 'C_alpha', 'C_carbo')
                             # include N, C_alpha, C_carbo in start and end.
        self.group2frag = dict(start = dict(N='y', C_alpha='z', C_carbo='x'),
                                 end = dict(N='c', C_alpha='a', C_carbo='b'))
        
        self.modifications = {(int(number) - 1, group): Formula(atom_cnt)
                              for number, mods in modifications.items()
                              for group, atom_cnt in mods.items()}
        
        self.formula = sum(self[number, group]
                           for number in range(len(self))
                           for group in self.groups)

        self.blocked_fragments = set(blocked_fragments)
        self.fragments         = fragments
        self.distance_charges  = int(distance_charges)

        if block_prolines:
            for i, f in enumerate(self.fasta):
                if f == 'P':
                    self.blocked_fragments.add('c' + str(i))
                    z_frag_No = len(self) - i
                    self.blocked_fragments.add('z' + str(z_frag_No))

    def _get_amino_acid(self, number, group):
        """Get amino acid of the precursor."""
        formula = self.amino_acids[(self.fasta[number], group)] +\
                  self.modifications.get((number, group), 0)
        # Modifying termini: Kaltashov O., Mass Spectrometry in Biophysics
        if number == 0 and group == 'N':
            formula['H'] += 1  #  H for N terminus
        if number == (len(self.fasta)-1) and group == 'C_carbo':
            formula['O'] += 1  # additional H and O
            formula['H'] += 1  # for the C terminus
        formula.check_positivity()
        return formula

    # TODO this should work with sequences, e.g. 1:10 
    def __getitem__(self, key):
        """Get amino acid of the precursor.

        Parameters
        ==========
        key : int or tuple(int, str)
            The string should be C_carbo, C_alpha, or N.
            The int describes the number of the amino acid, starting from
            zero, counting from N terminus to C terminus.
        Returns
        =======
        out : Formula
        """
        try:
            return self._get_amino_acid(*key)
        except TypeError:
            return sum(self._get_amino_acid(key, group)
                       for group in self.groups)
        except ValueError:
            raise KeyError("Supply '(number, group)' or just 'group'.")

    def __repr__(self):
        out = "({name} {q}+ {fasta}".format(**self.__dict__)
        out += '-modified)' if self.modifications else ')'
        return out

    def __len__(self):
        return len(self.fasta)

    def _protonate(self, frag):
        a, b, c = {'p': (1,  0, 1),
                   'c': (0, -1, 0),
                   'z': (0,  0, 1)}[frag]
        for q in range(1, self.q + a):
            for g in range(b, self.q - q + c):
                yield (q, g)

    def a_fragments(self):
        """Generate a fragments."""
        raise NotImplementedError

    def b_fragments(self):
        """Generate b fragments."""
        raise NotImplementedError

    def c_fragments(self):
        """Generate c fragments."""
        # 'H1' to be a 'c' fragment, not the H on the N terminus
        formula = Formula('H1')
        for number in range(len(self.fasta)):
            formula += self[number, 'N']
            name = 'c' + str(number)
            if name not in self.blocked_fragments:
                yield (name, formula.copy())
            formula += self[number, 'C_alpha']
            formula += self[number, 'C_carbo']

    def x_fragments(self):
        """Generate x fragments."""
        raise NotImplementedError

    def y_fragments(self):
        """Generate y fragments."""
        raise NotImplementedError

    def z_fragments(self):
        """Generate z fragments."""
        formula = Formula()
        for number in range(len(self)-1, -1, -1):
            formula += self[number, 'C_carbo']
            formula += self[number, 'C_alpha']
            side_chain_len = len(self) - number
            name = 'z' + str(side_chain_len)
            if name not in self.blocked_fragments:
                yield (name, formula.copy())
            formula += self[number, 'N']

    def uncharged_molecules(self):
        """Generate uncharged molecules."""
        yield ('precursor', self.formula.copy())
        for frag_type in self.fragments:
            frags = getattr(self, frag_type + '_fragments')()
            for frag in frags:
                yield frag

    def molecules(self):
        """Generate charged molecules."""
        for name, formula in self.uncharged_molecules():
            if name[0] != 'p':
                side_chain_len = int(name[1:])
            else:
                side_chain_len = len(self)
            #TODO: protonation rules other than for cz
            for q, g in self._protonate(name[0]):
                potential_charges_cnt = \
                    side_chain_len // self.distance_charges
                if side_chain_len % self.distance_charges > 0:
                    potential_charges_cnt += 1
                    # +0000 +0000 00+  at most 3 charges
                if potential_charges_cnt >= q:
                    mol = Molecule(formula  = formula,
                                   iso_calc = self.iso_calc,
                                   q        = q,
                                   g        = g)
                    yield mol, name

    def __hash__(self):
        """Get a hash from the precursor's unique id.

        Unique id consists of a name, fasta, charge, and modifications.
        """
        return hash((self.name,
                     self.fasta,
                     self.q,
                     flatten_modification(self.modifications)))

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            A = self.name == other.name
            B = self.fasta == other.fasta
            C = self.q == other.q
            return A and B and C
        else:
            return False


# todo: a fishy class. Have to do something about it.
class FalsePrecursor(Molecule):
    """This is a class needed only for the Molecules class.

    That's obiously fishy and we should inspect, if it is necessary at all.
    """
    def __init__(self,
                 name,
                 formula,
                 iso_calc,
                 q        = 0):
        self.name      = name
        self.formula   = Formula(formula)
        self.q         = int(q)
        self.intensity = 0.0
        self.iso_calc  = iso_calc
        self.real      = False

    def __hash__(self):
        return hash((self.name,
                     self.formula.str_with_charges(self.q),
                     self.q))

    def __eq__(self, other):
        A = self.name    == other.name
        B = self.formula == other.formula
        C = self.q       == other.q
        return A and B and C

    def __repr__(self):
        return "({name} q={q} I={I_int})".format(
            I_int = int(self.intensity),
            **self.__dict__)


def precursor(fasta,
              q,
              iso_calc,
              name              = "",
              modifications     = {},
              fragments         = "cz",
              blocked_fragments = set(['c0']),
              block_prolines    = True,
              distance_charges  = 5):
    """Prepare a ready precursor."""
    prec = Precursor(fasta,
                     q,
                     iso_calc,
                     name,
                     modifications,
                     fragments,
                     blocked_fragments,
                     block_prolines,
                     distance_charges)
    return prec