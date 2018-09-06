from   bisect import bisect_left
import csv
import numpy as np
from   operator import itemgetter

from masstodon.data.constants import infinity
from masstodon.misc.strings   import repr_long_list
from masstodon.parse.path     import parse_path
# to do: add bokeh as other possible package.

class Measure(object):
    """Store a discrete finite measure with atoms on the real line."""

    def __init__(self, atoms=np.array([]),
                       masses=np.array([]),
                       sort=True):
        """Initialize a measure.

        Parameters
        ==========
        atoms : numpy array
            The atoms upon which the measure holds the mass.
        masses : numpy array
            The masses on atoms.

        Remarks
        =======
        If you have sorted the measure beforehand, do set 'is_sorted' to True.

        """
        self.atoms  = atoms
        self.masses = masses
        if sort:
            self.sort()
        self._store_names = ('atom', 'mass')

    def sort(self):
        """Sort measure by atomic values."""
        atom_sorted = np.argsort(self.atoms)
        self.atoms = self.atoms[atom_sorted]
        self.masses = self.masses[atom_sorted]


    def __has_type_of(self, other):
        """Assert that 'self' and 'other' have the same type."""
        assert self.__class__.__name__ == other.__class__.__name__,\
             "\tIllegal to add class {0} to class {1}.\n".format(
                other.__class__.__name__,\
                self.__class__.__name__)

    def __add__(self, other):
        """Add two measures.

        Parameters
        ----------
        other : Measure
            A measure we want to stack on top of this one.

        """
        atoms = np.concatenate((self.atoms, other.atoms))
        masses = np.concatenate((self.masses, other.masses))
        self.__has_type_of(other)
        new_measure = self.__class__(atoms, masses)
        new_measure.__aggregate()
        return new_measure

    def __radd__(self, other):
        """Add two measures.

        Parameters
        ----------
        other : Measure
            A measure we want to stack on top of this one.

        """
        if other == 0:
            return self
        else:
            return self.__add__(other)

    def __iadd__(self, other):
        """Add two measures.

        Parameters
        ----------
        other : Measure
            A measure we want to stack on top of this one.

        """
        self.__has_type_of(other)
        self.atoms  = np.concatenate((self.atoms, other.atoms))
        self.masses = np.concatenate((self.masses, other.masses))
        self.__aggregate()
        return self

    def copy(self):
        """Make a deep copy of me."""
        out = self.__class__(self.atoms, self.masses)
        return out

    def __mul__(self, scalar):
        """Multiply by a scalar."""
        if scalar == 0:
            return self.__class__()
        elif scalar == 1:
            return self.copy()
        else:
            return self.__class__(self.atoms, scalar * self.masses)

    def __rmul__(self, scalar):
        """Multiply by a scalar."""
        if scalar == 1:
            return self
        else:
            return self.__mul__(scalar)

    def __imul__(self, scalar):
        """Multiply by a scalar."""
        if scalar != 1:
            self.masses = self.masses * scalar
        return self

    def __aggregate(self):
        """Aggregate masses with the same atoms."""
        self.atoms, indices = np.unique(self.atoms, return_inverse=True)
        self.masses = np.bincount(indices, weights=self.masses)

    def round_atoms(self, precision=infinity):
        """Round the atoms of the measure to a given precision.

        Parameters
        ----------
        precision : integer
            The number of digits after which the atoms' masses get rounded.
            E.g. if set to 2, then number 3.141592 will be rounded to 3.14.
            Defaults to 'inf', which prevents any rounding.

        """
        if precision != infinity:
            self.atoms = np.around(self.atoms, precision)
            self.__aggregate()

    def trim(self, cut_off):
        """Trim masses below the provided cut off.

        Parameters
        ----------
        cut_off : float

        """
        if cut_off > 0:
            self.atoms  = self.atoms[self.masses >= cut_off]
            self.masses = self.masses[self.masses >= cut_off]

    def split_measure(self, cut_off):
        """Split measure into two according to the cut off on masses.

        Retain the measure with masses greater or equal to the cut off.
        Parameters
        ----------
        cut_off : float
        Returns
        ----------
        other : Measure
            A measure with masses strictly below the cut off.
        """
        other = self.__class__(self.atoms[self.masses < cut_off],
                               self.masses[self.masses < cut_off])
        self.trim(cut_off)
        return other

    def get_P_set_cut_off(self, P=.99):
        """Get the cut off resulting in optimal P-set.

        Parameters
        ----------
        P : float
            The percentage of the initial total value that the new measure
            will contain. The new measure contains only atoms with
            heighest masses.

        """
        # TODO replace this with something brighter, as linearly operating heap.
        assert 0.0 <= P <= 1.0, "Wrong P for P-optimal set."
        total_value = self.masses.sum()
        i = 0
        S = 0.0
        masses = np.sort(self.masses)[::-1]
        for intensity in masses:
            S += intensity
            if S < P:
                break
        return intensity

    def __repr__(self):
        """Represent the measure."""
        if len(self.atoms) > 0:
            out = "{0}:\n\t{1} = {2}\n\t{3} = {4}\n".format(self.__class__.__name__,
                                                            self._store_names[0],
                                                            repr_long_list(self.atoms)[1:-1],
                                                            self._store_names[1],
                                                            repr_long_list(self.masses)[1:-1])
        else:
            out = self.__class__.__name__ + " (empty)"
        return out

    def __len__(self):
        """Get size of the measure: the number of atoms."""
        return len(self.atoms)

    def __iter__(self):
        """Iterate over pairs (atom, mass)."""
        return zip(self.atoms, self.masses)

    #TODO: 
    # extend to any slice
    def __getitem__(self, key):
        """Filter atoms between 'L' and 'R'.

        Parameters
        ==========
        key : tuple
            Either (left_mz, right_mz) or
            (left_mz, right_mz, provide_idx),
            where 'provide_idx' is boolean.
        Returns
        =======
        out : generator
            Generate tuples '(mz, intensity)'
            or '(idx, mz, intensity)',
            where 'idx' is the unique ID of the atom.

        """
        try:
            if len(key) == 2:
                L, R = key
                provide_idx = False
            elif len(key) == 3:
                L, R, provide_idx = key
            idx = bisect_left(self.atoms, L)
            while self.atoms[idx] <= R:
                if not provide_idx:
                    yield self.atoms[idx], self.masses[idx]
                else:
                    yield idx, self.atoms[idx], self.masses[idx]
                idx += 1
        except IndexError:
            # WTF???
            return


    def total_mass(self):
        return self.masses.sum()

    def write(self, path):
        """Write the spectrum to a csv or tsv file.

        Parameters
        ==========
        path : str
            A path to the file to write to.
        """
        file_path, file_name, file_ext = parse_path(path)
        delimiter = ',' if file_ext == '.csv' else '\t'
        with open(path, 'w') as csvfile:
            writer = csv.writer(csvfile, delimiter=delimiter)
            writer.writerow(self._store_names)
            for atom, mass in self:
                writer.writerow([atom, mass])
