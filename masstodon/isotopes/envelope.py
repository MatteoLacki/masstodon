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
import numpy as np
try:
    import matplotlib.pyplot as plt
except RuntimeError:
    pass

from masstodon.measure.measure import Measure

# Compared to Measure, the multiplication changes atoms, rather then masses.
class Envelope(Measure):
    """Store an isotopic distribution."""

    def __init__(self,
                 mz          = np.array([]),
                 probability = np.array([]),
                 sort        = True):
        """Initialize an isotopic distribution.

        Parameters
        ----------
        mz : numpy array
            Mass to charge ratios of the isotopic distribution.
        probability : numpy array
            Probabilities of the isotopic distribution.

        """
        super().__init__(mz, probability, sort)
        self._store_names = ('m/z', 'probability')

    @property
    def mz(self):
        """Get mass over charge ratios"""
        return self.atoms

    @mz.setter
    def mz(self, mz):
        """Set m/z ratios."""
        self.atoms = mz

    @property
    def probability(self):
        """Get probabilities."""
        return self.masses

    @probability.setter
    def probability(self, probability):
        """Set probabilities."""
        self.masses = probability

    def copy(self):
        """Make a deep copy of me."""
        out = self.__class__(self.mz, self.probability)
        return out

    def __mul__(self, scalar):
        """Multiply by a scalar."""
        if scalar == 0:
            return self.__class__()
        elif scalar == 1:
            return self.copy()
        else:
            return self.__class__(scalar * self.mz, self.probability)

    def __rmul__(self, scalar):
        """Multiply by a scalar."""
        if scalar == 1:
            return self
        else:
            return self.__mul__(scalar)

    def __imul__(self, scalar):
        """Multiply by a scalar."""
        if scalar != 1:
            self.mz = self.mz * scalar
        return self

    def __truediv__(self, scalar):
        return self.__div__(scalar)

    def __div__(self, scalar):
        inverse = 1.0 / scalar
        return self.__mul__(inverse)

    def __rdiv__(self, scalar):
        inverse = 1.0 / scalar
        return self.__rmul__(inverse)

    def __idiv__(self, scalar):
        inverse = 1.0 / scalar
        return self.__imul__(inverse)

    def add_mass_divide_by_charge(self, mass, q):
        """Add mass and divide by charge.

        This avoids double copying.
        (of first adding and then dividing)
        """
        return self.__class__((self.mz + mass)/q, self.probability)

    def head_mz(self):
        return self.mz[0]

    def tail_mz(self):
        return self.mz[-1]

    def max_peak(self):
        i = np.argmax(self.mz)
        return self.mz[i], self.probability[i]

    def round_mz(self, precision):
        """Round the mass to charge ratios to a given precision.

        Parameters
        ----------
        precision : integer
            The number of digits after which the atoms' masses get rounded.
            E.g. if set to 2, then number 3.141592 will be rounded to 3.14.

        """
        self.round_atoms(precision)

    def plot(self,
             plt_style = 'dark_background',
             show      = True):
        plt.style.use(plt_style)
        plt.vlines(self.mz, [0], self.probability, colors='blue')
        if show:
            plt.show()


def envelope(mz, probability, sort):
    return Envelope(mz, probability, sort)