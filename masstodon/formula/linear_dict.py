from collections import MutableMapping

class LinearDict(MutableMapping):
    """A dictionary extended to linear operations."""
    def __init__(self, *args, **kwds):
        self._storage = dict(*args, **kwds)

    def __getitem__(self, key):
        return self._storage.get(key, 0)

    def __setitem__(self, key, value):
        self._storage[key] = value

    def __delitem__(self, key):
        del self._storage[key]

    def __iter__(self):
        return iter(self._storage)

    def __len__(self):
        return len(self._storage)

    def __add__(self, other):
        out = self.copy()
        out.__iadd__(other)
        return out

    def __radd__(self, other):
        if other == 0:
            return self
        else:
            return self.__add__(other)

    def __iadd__(self, other):
        if other != 0:
            try:
                for key in set().union(self._storage, other._storage):
                    self._storage[key] = self._storage.get(key, 0) + \
                                          other._storage.get(key, 0)
            except (TypeError, AttributeError):
                raise TypeError("Adding supported for dict-like structures.")
        return self

    def __repr__(self):
        return self._storage.__repr__()

    def copy(self):
        out = self.__class__(self)
        return out

    def __mul__(self, scalar):
        """Implement."""
        raise NotImplementedError