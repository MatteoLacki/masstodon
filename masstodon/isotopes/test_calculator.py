"""Testing the isotope calculator."""
from __future__ import absolute_import, division, print_function
import numpy as np
import unittest

from MassTodonPy.IsotopeCalculator.IsotopeCalculator import IsotopeCalculator
from MassTodonPy.Misc.binomial import binomial

class TestPeakPicker(unittest.TestCase):
    def test_isotope_calculator(self):
        print("Testing the get_envelope function.")

        iso_calc = IsotopeCalculator(mz_digits=0,
                                     _masses={'T': [1,1000]},
                                     _probabilities={'T': [0.5, 0.5]})

        distribution = iso_calc.get_envelope(formula='T10',
                                             joint_probability=1.0)
        # E_ = expected
        E_atoms = np.array([1.*(10-i) + 1000.*i for i in range(11)])
        self.assertTrue(all(E_atoms == distribution.atoms))
        # bin(n, k) p**k q**(N-k) = bin(n, k) 0.5**10
        E_masses = np.array([binomial(10,i) for i in range(11)]) / 2**10
        masses_comparison = list(np.abs(E_masses - distribution.masses))
        for mass_diff in masses_comparison:
            self.assertTrue(mass_diff < 1e-12)


if __name__ == "__main__":
    unittest.main()
