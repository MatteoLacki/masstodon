# TODO: this could be further optimized by the use of numpy
# 
# observation: the charges must be an input here,
# otherwise, we will not have an easy access to H mass.
from math import sqrt
import numpy as np

from masstodon.data.constants    import infinity
from masstodon.formula.formula   import as_formula
from masstodon.data.isotopes     import get_isotopic_masses, get_isotopic_probabilities
from masstodon.isotopes.misc     import get_mean_and_variance, cdata2numpyarray  # TODO IsoSpec 2.0
from masstodon.isotopes.envelope import envelope
from masstodon.isotopes.isospec  import isospec_numpy



class IsotopeCalculator(object):
    def __init__(self,
                 digits         = infinity,
                 _masses        = get_isotopic_masses(),
                 _probabilities = get_isotopic_probabilities(),
                 _isotope_DB    = {}):
        """Initialize the isotopic calculator."""
        self.digits         = digits
        self._masses        = _masses
        self._probabilities = _probabilities
        self._isotope_DB    = _isotope_DB
        self._mean_mass      = {}
        self._mean_variance  = {}
        self._lightiest_mass = {}
        self._heaviest_mass  = {}
        self._probabliest_mass = {} # no respect for grammar!!! ever again!!! 
        for e in self._probabilities:
            self._mean_mass[e], self._mean_variance[e] = \
                get_mean_and_variance(self._masses[e],
                                      self._probabilities[e])
            self._lightiest_mass[e] = min(self._masses[e])
            self._heaviest_mass[e]  = max(self._masses[e])
            self._probabliest_mass[e] = self._masses[e][np.argmax(self._probabilities[e])]

    def mass_to_mz(self, mass, q=0, g=0):
        """Adjust mass to m/z."""
        if q > 0:
            H_mass = self._masses['H'][0]
            mz = (mass + (q + g) * H_mass) / q
            return mz
        else:
            return mass

    def lightiest_mz(self, formula, q=0, g=0):
        mass = sum(self._lightiest_mass[e] * formula[e] for e in formula)        
        return self.mass_to_mz(mass, q, g)

    def heaviest_mz(self, formula, q=0, g=0):
        mass = sum(self._heaviest_mass[e] * formula[e] for e in formula)
        return self.mass_to_mz(mass, q, g)

    def most_probable_mz(self, formula, q=0, g=0):
        mass = sum(self._probabliest_mass[e] * formula[e] for e in formula)
        return self.mass_to_mz(mass, q, g)

    def monoisotopic_mz(self, formula, q=0, g=0):
        """Calculate monoisotopic m/z or mass of a molecule."""
        return self.lightiest_mz(formula, q, g)

    def mean_mz(self, formula, q=0, g=0):
        """Calculate average m/z or mass of a molecule."""
        mass = sum(self._mean_mass[e] * formula[e] for e in formula)
        return self.mass_to_mz(mass, q, g)

    def sd_mz(self, formula, q, g=0):
        """Calculate standard deviation of m/z or mass for a molecule."""
        var = sum(self._mean_variance[e] * formula[e] for e in formula)
        if q > 0:
            H_var = self._mean_variance['H']
            var += H_var * (q + g)
            return sqrt(var) / q
        else:
            return sqrt(var)

    def _make_envelope(self, formula, prob, sort = True):
        el_cnt  = tuple(formula.values())
        el_mass = tuple(self._masses[e] for e in formula)
        el_prob = tuple(self._probabilities[e] for e in formula)
        m, p    = isospec_numpy(el_cnt, el_mass, el_prob, prob)
        env     = envelope(m, p, sort=True)
        return env

    # add possibility to call IsoSpec each time for precise calculations.
    def __call__(self,
                 formula,
                 prob     = .999,
                 q        = 0,
                 g        = 0,
                _memoize  = False):
        """Get an isotopic envelope."""
        if _memoize:
            env_key = (str(formula), prob)
            try:
                env = self._isotope_DB[env_key]
            except KeyError:
                env = self._make_envelope(formula, prob)
                if self.digits < infinity:
                    # round to prescibed level.
                    env.round_mz(self.digits)
                self._isotope_DB[env_key] = env
        else:
            formula = as_formula(formula)
            env = self._make_envelope(formula, prob)
            if self.digits < infinity:
                # round to prescibed level.
                env.round_mz(self.digits)
        # simplification: the q and g are only shifting the
        # distribution, rather than affecting the whole distribution.
        if q > 0:
            H_mass = self._masses['H'][0]
            env = env.add_mass_divide_by_charge(H_mass * (q + g), q)

        env.round_mz(self.digits)
        return env


def isotope_calculator(digits         = infinity,
                       _masses        = get_isotopic_masses(),
                       _probabilities = get_isotopic_probabilities(),
                       _isotope_DB    = {}):
    IC = IsotopeCalculator(digits         = digits,
                           _masses        = _masses,
                           _probabilities = _probabilities,
                           _isotope_DB    = _isotope_DB)
    return IC
