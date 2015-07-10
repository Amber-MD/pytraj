from __future__ import absolute_import
import operator
from pytraj.utils import _import_numpy
from pytraj._cyutils import _fast_count

np = _import_numpy()[-1]


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

    """place holder for all cpptraj' output
    """

    def __init__(self, dset=None, copy=True):
        """
        Parameters
        ----------
        dset : Cpptraj's DataSet or a Dict

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
            self.legend = getattr(dset, 'legend', "")
            self.name = getattr(dset, 'name', "")
            self.aspect = getattr(dset, 'aspect', 'unknown')
            self.idx = getattr(dset, 'idx', 0)
            self.format = getattr(dset, 'format', None)
            self.scalar_type = getattr(dset, 'scalar_type', 'unknown')
            if hasattr(dset, 'cpptraj_dtype'):
                self.cpptraj_dtype = dset.cpptraj_dtype
            else:
                from pytraj.datasets import DataSet
                if isinstance(dset, DataSet):
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
        import numpy as np
        self._values = np.asarray(new_values)

    def __iter__(self):
        for x in self._values:
            yield x

    def __getitem__(self, idx):
        return self._values[idx]

    def __setitem__(self, idx, value):
        self._values[idx] = value

    @property
    def size(self):
        return len(self._values)

    @property
    def dtype(self):
        return self._values.dtype

    def astype(self, t):
        self._values = self._values.astype(t)

    @property
    def key(self):
        return self.legend

    @key.setter
    def key(self, new_key):
        self.legend = new_key

    @property
    def data(self):
        return self.values

    def __str__(self):
        size = self.size
        key = self.key
        dtype = self.dtype
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

    def shallow_copy(self):
        """everything is copied but `self.values`
        """
        return self.__class__(self, copy=False)

    def is_empty(self):
        return len(self.values) == 0

    def append(self, value, axis=None):
        self.values = np.append(self.values[:], value, axis=axis)

    @property
    def ndim(self):
        return self.values.ndim

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

    def to_pyarray(self):
        from array import array
        if 'int' in self.dtype.name:
            return array('i', self.values.flatten())
        else:
            return array('d', self.values.flatten())

    def count(self, value=None):
        """
        Parameters
        value : int, optional

        Examples
        --------
        ds.count()
        ds.count(1)
        """
        if value is None:
            from collections import Counter
            return Counter(self.values)
        else:
            return _fast_count(self.values, value)

    def hist(self, plot=True, show=True, *args, **kwd):
        """
        Parameters
        ----------
        plot : bool, default True
            if True, use `matplotlib` to plot. 
            if False, return `2D numpy array`
        """
        if not plot:
            import numpy as np
            return np.histogram(self.values)
        else:
            try:
                from matplotlib import pyplot as plt
                ax = plt.hist(self, *args, **kwd)
                if show:
                    plt.show()
                return ax
            except ImportError:
                raise ImportError("require matplotlib")

    def split(self, n_chunks_or_array):
        """split `self.data` to n_chunks

        Notes : require numpy (same as `array_split`)
        """
        return np.array_split(self.to_ndarray(), n_chunks_or_array)

    def plot(self, show=True, *args, **kwd):
        """return matplotlib object
        Notes
        ----
        Need to over-write this method for subclass if needed.
        """
        from pytraj.utils import _import
        try:
            from matplotlib import pyplot as plt
            ax = plt.subplot(111)
            if self.ndim == 2 and 'B-factor' in self.key:
                values = self.values.T
                # good for bfactor plotting
                ax.plot(values[0], values[1], *args, **kwd)
            else:
                ax.plot(self, *args, **kwd)
            if show:
                plt.show()
            return ax
        except ImportError:
            raise ImportError("require matplotlib")

    def chunk_average(self, n_chunk):
        from pytraj.tools import chunk_average
        return chunk_average(self, n_chunk)

    def std(self, *args, **kwd):
        import numpy as np
        return np.std(self.values, *args, **kwd)

    def topk(self, k):
        """pick top k max-values
        Returns
        -------
        a list with len = k

        # TODO : array?
        """
        return sorted(self.values, reverse=True)[:k]

    def head(self, k=20, restype='ndarray'):
        if restype == 'ndarray':
            return self.values[:k]
        elif restype == 'list':
            return self.tolist()[:k]

    def tail(self, k=20):
        return self.values[-k:]

    def filter(self, func):
        """return a numpy array with all elements that satisfy `func`

        Example
        -------
        >>> d0 = traj.calc_radgyr(dtype='dataset')[0]
        >>> d0.filter(lambda x : 105. < x < 200.)
        """
        import numpy as np
        return np.array(list(filter(func, self.values)))

    def _normalize(self):
        """simple normalization
        (data - data.mean()) / (data.max() - data.min())
        """
        return (self - self.mean()) / (self.max() - self.min())

    def sum(self, axis=None, *args, **kwd):
        return self.values.sum(axis=axis, *args, **kwd)

    def avg(self, *args, **kwd):
        # use *args, **kwd so we can use numpy too
        return self.mean()

    def mean(self, *args, **kwd):
        # use *args, **kwd so we can use numpy too
        return self.values.mean()

    def min(self, *args, **kwd):
        # use *args, **kwd so we can use numpy too
        return self.values.min()

    def max(self, *args, **kwd):
        # use *args, **kwd so we can use numpy too
        return self.values.max()

    def median(self, *args, **kwd):
        # use *args, **kwd so we can use numpy too
        return np.median(self.values)

    def flatten(self):
        return self.values.flatten()

    # below are copied from `dask` package: New BSD
    # see pytraj/licenses/externals/dask.txt for license
    def __abs__(self):
        return elemwise(operator.abs, self)

    def __add__(self, other):
        return elemwise(operator.add, self, other)

    def __radd__(self, other):
        return elemwise(operator.add, other, self)

    def __and__(self, other):
        return elemwise(operator.and_, self, other)

    def __rand__(self, other):
        return elemwise(operator.and_, other, self)

    def __div__(self, other):
        return elemwise(operator.div, self, other)

    def __rdiv__(self, other):
        return elemwise(operator.div, other, self)

    def __invert__(self):
        return elemwise(operator.invert, self)

    def __lshift__(self, other):
        return elemwise(operator.lshift, self, other)

    def __rlshift__(self, other):
        return elemwise(operator.lshift, other, self)

    def __mod__(self, other):
        return elemwise(operator.mod, self, other)

    def __rmod__(self, other):
        return elemwise(operator.mod, other, self)

    def __mul__(self, other):
        return elemwise(operator.mul, self, other)

    def __rmul__(self, other):
        return elemwise(operator.mul, other, self)

    def __neg__(self):
        return elemwise(operator.neg, self)

    def __or__(self, other):
        return elemwise(operator.or_, self, other)

    def __pos__(self):
        return self

    def __ror__(self, other):
        return elemwise(operator.or_, other, self)

    def __rpow__(self, other):
        return elemwise(operator.pow, other, self)

    def __rshift__(self, other):
        return elemwise(operator.rshift, self, other)

    def __rrshift__(self, other):
        return elemwise(operator.rshift, other, self)

    def __sub__(self, other):
        return elemwise(operator.sub, self, other)

    def __rsub__(self, other):
        return elemwise(operator.sub, other, self)

    def __truediv__(self, other):
        return elemwise(operator.truediv, self, other)

    def __rtruediv__(self, other):
        return elemwise(operator.truediv, other, self)

    def __floordiv__(self, other):
        return elemwise(operator.floordiv, self, other)

    def __rfloordiv__(self, other):
        return elemwise(operator.floordiv, other, self)

    def __xor__(self, other):
        return elemwise(operator.xor, self, other)

    def __rxor__(self, other):
        return elemwise(operator.xor, other, self)
    # end of copy from dask
