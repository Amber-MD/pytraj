from __future__ import absolute_import
import operator
import numpy as np
from pytraj._cyutils import _fast_count
from pytraj.datasets import Dataset


def elemwise(op, self, other=None):
    if other:
        if hasattr(other, 'values'):
            _other = other.values
        else:
            _other = other
        if hasattr(self, 'values'):
            _self = self.values
        else:
            _self = self
        return op(_self, _other)
    else:
        return op(self.values)


class DataArray(object):
    """place holder for all cpptraj' output.

    Examples
    --------
    >>> from pytraj.array import DataArray
    >>> arr0 = DataArray([0, 2, 4])
    >>> arr1 = DataArray({'x' : [0, 2, 4]})
    >>> print(arr1)
    <pytraj.array.DataArray: size=3, key=x, dtype=int64, ndim=1> 
    values:
    [0 2 4]
    >>> print(arr1.values)
    """

    def __init__(self, dset=None, copy=True):
        """
        Parameters
        ----------
        dset : Cpptraj's Dataset or a Dict

        Examples
        --------
        >>> DataArray({'x' : [3, 5, 6]])
        """
        if isinstance(dset, dict):
            assert len(dset.keys()) == 1, "single dict"
            key = list(dset.keys())[0]
            self.key = key
            if copy:
                self._values = np.asarray(dset[key]).copy()
            else:
                self._values = np.asarray(dset[key])
            self.name = ""
            self.aspect = "unknown"
            self.idx = 0
            self.format = None
            self.scalar_type = 'unknown'
            self.cpptraj_dtype = None
        else:
            self.key= getattr(dset, 'key', "")
            self.name = getattr(dset, 'name', "")
            self.aspect = getattr(dset, 'aspect', 'unknown')
            self.idx = getattr(dset, 'idx', 0)
            self.format = getattr(dset, 'format', None)
            self.scalar_type = getattr(dset, 'scalar_type', 'unknown')
            if hasattr(dset, 'cpptraj_dtype'):
                self.cpptraj_dtype = dset.cpptraj_dtype
            else:
                if isinstance(dset, Dataset):
                    self.cpptraj_dtype = dset.dtype
                else:
                    self.cpptraj_dtype = None

        if dset is not None and not isinstance(dset, dict):
            if hasattr(dset, 'values'):
                values = np.asarray(dset.values)
            else:
                values = np.asarray(dset)
            if copy:
                self._values = values.copy()
            else:
                self._values = values

    @classmethod
    def from_dict(cls, d):
        assert isinstance(d, dict), "must be a dict"
        return cls(d)

    @property
    def values(self):
        return self._values

    @values.setter
    def values(self, new_values):
        self._values = np.asarray(new_values)

    @property
    def ndim(self):
        return self.values.ndim

    def __iter__(self):
        for x in self._values:
            yield x

    def __getitem__(self, idx):
        return self._values[idx]

    def __setitem__(self, idx, value):
        self._values[idx] = value

    def transpose(self):
        d = self.__class__(self, copy=False)
        d.values = d.values.T
        return d

    T = property(transpose)

    @property
    def size(self):
        return len(self._values)

    @property
    def dtype(self):
        return self._values.dtype

    def astype(self, t):
        self._values = self._values.astype(t)

    @property
    def data(self):
        return self.values

    def __str__(self):
        size = self.size
        key = self.key
        msg0 = """<pytraj.array.DataArray: size={0}, key={1}, dtype={2}, ndim={3}> """.format(
            size, key, self.dtype, self.ndim)
        value_str = self.values.__str__()
        return msg0 + '\nvalues:\n' + value_str

    def __repr__(self):
        return self.__str__()

    def __array__(self):
        return self.values

    def __len__(self):
        return len(self.values)

    def copy(self):
        """deep copy"""
        new_ds = self.__class__(self)
        new_ds.values = self.values.copy()
        return new_ds

    def _shallow_copy(self):
        """everything is copied but `self.values`
        """
        return self.__class__(self, copy=False)

    def is_empty(self):
        return len(self.values) == 0

    def append(self, value, axis=None):
        self.values = np.append(self.values[:], value, axis=axis)

    @property
    def shape(self):
        return self.values.shape

    def tolist(self):
        return self.values.tolist()

    def to_ndarray(self, copy=False):
        if copy:
            return self.values.copy()
        else:
            return self.values

    def to_dict(self):
        return {self.key: self.values}

    def flatten(self):
        return self.values.flatten()
