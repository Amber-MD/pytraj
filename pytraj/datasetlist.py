from __future__ import absolute_import
from array import array
import numpy as np
from functools import partial
from pytraj.datasets.DatasetList import DatasetList as DSL
from pytraj.externals._json import to_json, read_json
from pytraj.externals._pickle import to_pickle, read_pickle
from pytraj.utils import _import_pandas, is_int, is_array, is_generator
from pytraj.compat import string_types, callable
from pytraj.datafiles import DataFile
from pytraj.core.cpp_core import ArgList
from pytraj.compat import map, iteritems
from pytraj.array import DataArray

__all__ = ['load_datafile', 'stack', 'DatasetList', 'from_pickle', 'from_json']


def _groupby(self, key):
    # adapted from `toolz` package.
    # see license in $PYTRAJHOME/licenses/externals/toolz.txt
    import collections
    if not callable(key):
        key = getter(key)
    d = collections.defaultdict(lambda: self.__class__().append)
    for item in self:
        d[key(item)](item)
    rv = {}
    for k, v in iteritems(d):
        rv[k] = v.__self__
    return rv


def from_pickle(filename):
    dslist = DatasetList()
    dslist.from_pickle(filename)
    return dslist


def from_json(filename):
    dslist = DatasetList()
    dslist.from_json(filename)
    return dslist


def from_dict(d):
    return DatasetList(d)


def load_datafile(filename):
    """load cpptraj's output"""
    ds = DatasetList()
    ds.read_data(filename)
    return ds


def _from_full_dict(full_dict):
    return DatasetList()._from_full_dict(full_dict)


def from_sequence(seq, copy=True):
    return DatasetList().from_sequence(seq, copy=copy)


def stack(args):
    """return a new DatasetList by joining (vstack)

    Parameters
    ----------
    args : list/tuple of DatasetList

    Notes
    -----
        similiar to numpy.vstack

    Examples
    --------
        d1 = calc_dssp(traj1, dtype='dataset')
        d2 = calc_dssp(traj2, dtype='dataset')
        d3 = stack((d1, d2))
    """
    is_subcriptable = not (isinstance(args, map) or is_generator(args))

    if not isinstance(args, (list, tuple, map)) and not is_generator(args):
        raise ValueError("must a tuple/list/map/generator")

    if is_subcriptable:
        dslist0 = args[0].copy()
    else:
        dslist0 = next(args)

    dslist_iter = args[1:] if is_subcriptable else args

    for dslist in dslist_iter:
        for d0, d in zip(dslist0, dslist):
            if d0.dtype != d.dtype:
                raise TypeError("Dont support stack different dtype together")
            if d0.key != d.key:
                raise KeyError("Don't support stack different key")
            d0.append(d.copy())
    return dslist0


concat_datasetlist = stack

from collections import OrderedDict


class _OrderedDict(OrderedDict):
    @property
    def values(self):
        from pytraj.tools import dict_to_ndarray
        arr = dict_to_ndarray(self)
        return np.array([[key for key in self.keys()], arr], dtype='object')

    def to_ndarray(self, with_key=False):
        from pytraj.tools import dict_to_ndarray
        arr = dict_to_ndarray(self)

        if not with_key:
            return arr
        else:
            return np.array([[key for key in self.keys()], arr])


class DatasetList(list):
    '''similiar to python's list but the data is labeled.
    Think as a dict-like and list-like object. This class is suitable for small
    datasets. For high performance, user should use pandas' DataFrame.

    Examples
    --------
    >>> import pytraj as pt
    >>> dslist = pt.multidihedral(traj, 'phi psi', dtype='dataset')
    >>> print(dslist)
    >>> print(dslist[0])
    >>> print(dslist[::2])
    >>> print(dslist[1::2])
    >>> print(dslist.values) # return raw numpy array
    >>> print(dslist.keys)
    >>> # dict-like
    >>> dslist['phi_2']
    '''

    def __init__(self, dslist=None, copy=False):
        if dslist:
            if isinstance(dslist, dict):
                # {'x': [1, 3, 5], 'y': [4, 7, 8]}
                for key, values in iteritems(dslist):
                    self.append(DataArray({key: values}), copy=copy)
            else:
                for d0 in dslist:
                    # always make a copy
                    # from DataArray(d0)
                    # so we set copy=False here
                    # to avoid copying twice
                    self.append(DataArray(d0), copy=copy)

    def copy(self):
        dslist = self.__class__()
        for d0 in self:
            dslist.append(d0, copy=True)
        return dslist

    def from_pickle(self, filename):
        ddict = read_pickle(filename)
        self._from_full_dict(ddict)

    def from_json(self, filename):
        ddict = read_json(filename)
        self._from_full_dict(ddict)

    def to_pickle(self, filename, use_numpy=True):
        to_pickle(self._to_full_dict(use_numpy), filename)

    def to_json(self, filename, use_numpy=True):
        full_dict = self._to_full_dict(use_numpy=use_numpy)
        for key in self.keys():
            d = full_dict[key]['values']
            if hasattr(d, 'dtype') and 'int' in d.dtype.name:
                full_dict[key]['values'] = d.tolist()
        to_json(full_dict, filename)

    def _from_full_dict(self, ddict):
        from pytraj.array import DataArray
        da = DataArray()

        if not isinstance(ddict, dict):
            raise ValueError("must be a dict")
        try:
            ordered_keys = ddict['ordered_keys']
        except KeyError:
            ordered_keys = ddict.keys()

        for key in ordered_keys:
            d = ddict[key]
            da.values = np.array(d['values'])
            da.aspect = d['aspect']
            da.name = d['name']
            da.idx = d['idx']
            da.key = key
            da.cpptraj_dtype = d['cpptraj_dtype']
            self.append(da)
        return self

    def _to_full_dict(self, use_numpy=True):
        """
        """
        ddict = {}
        ddict['ordered_keys'] = []
        for d in self:
            ddict['ordered_keys'].append(d.key)
            ddict[d.key] = {}
            _d = ddict[d.key]
            if use_numpy:
                _d['values'] = d.values
            else:
                _d['values'] = list(d.values)
            _d['name'] = d.name
            _d['cpptraj_dtype'] = d.cpptraj_dtype
            _d['aspect'] = d.aspect
            _d['idx'] = d.idx
        return ddict

    def hist(self, plot=True, show=True, *args, **kwd):
        """
        Parameters
        ----------
        plot : bool, default False
            if False, return a dictionary of 2D numpy array
            if True, return a dictionary of matplotlib object
        """
        hist_dict = dict(map(
            lambda x: (x.key, x.hist(plot=plot, show=False, *args, **kwd)),
            self))

        if show:
            # only show once
            from pytraj.plot import show
            show()
        return hist_dict

    def dtypes(self):
        return self.get_dtypes()

    def aspects(self):
        return self.get_aspects()

    def pipe(self, *funcs):
        """apply a series of functions to self's data
        """
        values = self.values
        for func in funcs:
            values = func(values)
        return values

    def apply(self, func):
        """update self's values from `funcs` and return `self`
        """
        for d0 in self:
            func(d0)
        return self

    def to_pyarray(self):
        if self.size > 1:
            raise NotImplementedError("only use `to_pyarray` for DataSet_1D")

        return self[0].to_pyarray()

    def __str__(self):
        safe_msg = "<pytraj.DatasetList with %s datasets>\n" % self.size
        if self.size == 0:
            return safe_msg
        msg = "\n\n".join("\n".join((d.key, d.values.__str__())) for d in self)
        str_first_3 = "\n\n".join("\n".join((d.key, d.values.__str__()))
                                  for d in self[:3])
        str_last_2 = "\n\n".join("\n".join((d.key, d.values.__str__()))
                                 for d in self[-2:])

        if self.size <= 5:
            return safe_msg + msg
        else:
            return safe_msg + str_first_3 + "\n...\n\n" + str_last_2

    def __repr__(self):
        return self.__str__()

    def __call__(self, *args, **kwd):
        return self.filter(*args, **kwd)

    def clear(self):
        self = []

    def is_empty(self):
        return self != []

    @property
    def size(self):
        return len(self)

    def __getslice__(self, i, j):
        return self.__getitem__(slice(i, j))

    def __getitem__(self, idx):
        """return a DataSet instance
        Memory view is applied (which mean this new insance is just alias of self[idx])
        Should we use a copy instead?
        """
        if self.size == 0:
            raise ValueError("size = 0: can not index")

        if is_int(idx):
            return super(DatasetList, self).__getitem__(idx)
        elif isinstance(idx, string_types):
            for d0 in self:
                if d0.key.upper() == idx.upper():
                    d0._base = self
                    return d0
        elif isinstance(idx, slice):
            # return new view of `self`
            start, stop, step = idx.indices(self.size)
            new_dslist = self.__class__()
            for _idx in range(start, stop, step):
                new_dslist.append(self[_idx], copy=False)
            return new_dslist
        elif is_array(idx) or isinstance(idx, list) and not isinstance(
                idx[0], bool):
            new_dslist = self.__class__()
            for _idx in idx:
                new_dslist.append(self[_idx], copy=False)
            return new_dslist
        elif isinstance(idx, tuple) and len(idx) == 2:
            return self[idx[0]][idx[1]]
        else:
            raise ValueError()

    def keys(self):
        """return a list"""
        tmp_list = []
        for d0 in self:
            tmp_list.append(d0.key)
        return tmp_list

    def get_aspects(self, is_set=True):
        """return a set of uniqure aspects if "is_set" = True
        else: return a full list
        """

        tmp_list = []
        for d0 in self:
            tmp_list.append(d0.aspect)
        if is_set:
            return set(tmp_list)
        else:
            return tmp_list

    def get_scalar_types(self):
        """return a list"""
        tmp_list = []
        for d0 in self:
            tmp_list.append(d0.scalar_type)
        return tmp_list

    def get_scalar_modes(self):
        """return a list"""
        tmp_list = []
        for d0 in self:
            tmp_list.append(d0.scalar_mode)
        return tmp_list

    def get_dtypes(self):
        """return a list"""
        tmp_list = []
        for d0 in self:
            tmp_list.append(d0.dtype)
        return tmp_list

    def iteritems(self):
        for key in self.keys():
            yield key, self[key]

    def items(self):
        return self.iteritems()

    def map(self, func):
        for d0 in self:
            yield func(d0)

    def filter(self, func, copy=False, *args, **kwd):
        """return a new view of DatasetList of func return True"""
        dslist = self.__class__()

        if isinstance(func, (string_types, list, tuple)):
            return self.grep(func, *args, **kwd)
        elif callable(func):
            for d0 in self:
                if func(d0, *args, **kwd):
                    dslist.append(d0, copy=copy)
            return dslist
        else:
            raise NotImplementedError("func must be a string or callable")

    def grep(self, key, mode='key', copy=False):
        """"return a new DatasetList object as a view of `self`

        Parameters
        ----------
        key : str or list
            keyword for searching
        mode: str, default='key'
            mode = 'key' | 'name' | 'dtype' | 'aspect'
        """
        import re

        # use __class__ so we can `filter` return the same class
        # we subclass this Cython class to python level
        dtmp = self.__class__()

        # dont free mem here
        for d0 in self:
            att = getattr(d0, mode)
            if isinstance(key, string_types):
                if re.search(key, att):
                    dtmp.append(d0, copy=copy)
            elif isinstance(key, (list, tuple)):
                for _key in key:
                    if re.search(_key, att):
                        dtmp.append(d0, copy=copy)
            else:
                raise ValueError("support string or list/tuple of strings")
        return dtmp

    def tolist(self):
        """return a list of list/array"""
        try:
            return [d0.tolist() for d0 in self]
        except:
            raise NotImplementedError("dont know how to convert to list")

    def to_dict(self, use_numpy=True, ordered_dict=False):
        """return a dict object with key=key, value=list"""
        _dict = dict
        if ordered_dict:
            # use OrderedDict
            _dict = _OrderedDict
        if use_numpy:
            return _dict((d0.key, d0.to_ndarray()) for d0 in self)
        else:
            return _dict((d0.key, d0.tolist()) for d0 in self)

    @property
    def values(self):
        """return read-only ndarray"""
        try:
            return self.to_ndarray()
        except:
            raise ValueError("don't know how to cast to numpy array")

    def to_ndarray(self, with_key=False):
        """
        Notes: require numpy
        """
        try:
            if self.size == 1:
                arr = self[0].values
            else:
                # more than one set
                arr = np.asarray([d0.values for d0 in self])
            if not with_key:
                return arr
            else:
                key_arr = np.array([d.key for d in self])
                return np.array([key_arr, arr.T])
        except:
            raise ValueError("don't know how to convert to ndarray")

    def to_dataframe(self):
        """return pandas' DataFrame

        Requires
        --------
        pandas
        """
        import pandas
        dict = _OrderedDict
        my_dict = dict((d0.key, d0.to_ndarray(copy=True)) for d0 in self)
        return pandas.DataFrame(my_dict)

    def mean(self, *args, **kwd):
        dict = _OrderedDict
        return dict((x.key, x.mean()) for x in self).values

    def median(self, *args, **kwd):
        """
        Notes: require numpy
        """
        dict = _OrderedDict
        return dict((x.key, x.median()) for x in self).values

    def std(self, axis=1):
        """
        Notes: require numpy
        """
        dict = _OrderedDict
        return dict((x.key, x.std()) for x in self).values

    def min(self, *args, **kwd):
        dict = _OrderedDict
        return dict((x.key, x.min()) for x in self).values

    def max(self, *args, **kwd):
        dict = _OrderedDict
        return dict((x.key, x.max()) for x in self).values

    def sum(self, restype='dict'):
        """
        Notes: require numpy
        """
        dict = _OrderedDict
        if restype == 'dict':
            return dict((x.key, x.sum()) for x in self).values
        elif restype == 'ndarray':
            return np.array([self.keys(), dict((x.key, x.sum())
                                               for x in self).to_ndarray()]).T

    def cumsum(self, axis=1):
        """Return the cumulative sum of the elements along a given axis.
        (from numpy doc)
        """
        dict = _OrderedDict
        return dict((x.key, np.cumsum(x.values)) for x in self).values

    def mean_with_error(self, other):
        ddict = defaultdict(tuple)
        for key, dset in self.iteritems():
            ddict[key] = dset.mean_with_error(other[key])
        return ddict

    def count(self, number=None):
        dict = _OrderedDict
        return dict((d0.key, d0.count(number)) for d0 in self)

    def read_data(self, filename, arg=""):
        df = DataFile()
        dslist = DSL()
        df.read_data(filename, ArgList(arg), dslist)
        df.from_sequence(dslist, copy=False)

    # pandas related
    def describe(self):
        import pandas as pd
        return self.to_dataframe().describe()

    def savetxt(self, filename='dslist_default_name.txt', labels=None):
        """just like `numpy.savetxt`
        Notes: require numpy
        """
        if labels is None:
            headers = "\t".join([d.key for d in self])
            headers = "frame\t" + headers
        else:
            headers = "frame\t" + labels

        frame_number = np.arange(self[0].size)
        # transpose `values` first
        values = np.column_stack((frame_number, self.values.T))
        formats = ['%8i'] + [d.format for d in self]
        np.savetxt(filename, values, fmt=formats, header=headers)

    def plot(self,
             show=True,
             use_seaborn=True,
             x=None,
             y=None,
             keys=[],
             autoset=True,
             xlim=None,
             ylim=None, *args, **kwd):
        """very simple plot for quickly visualize the data

        >>> dslist[['psi:7', 'phi:7']].plot()
        >>> dslist[['psi:7', 'phi:7']].plot(show=True)
        """
        if autoset:
            keys = self.keys() if not keys else keys

        if use_seaborn:
            try:
                import seaborn as snb
                snb.set()
            except ImportError:
                print("no seaborn package. skip importing")
        try:
            from matplotlib import pyplot as plt
            ax = plt.subplot(111)
            for idx, d0 in enumerate(self):
                lb = keys[idx] if keys else None
                if lb:
                    if d0.ndim == 2:
                        ax.plot(d0[0], d0[1], label=lb, *args, **kwd)
                    else:
                        ax.plot(d0, label=lb, *args, **kwd)
                else:
                    if d0.ndim == 2 and 'B-factors' in d0.key:
                        values = d0.values.T
                        ax.plot(values[0], values[1], *args, **kwd)
                    else:
                        ax.plot(d0, *args, **kwd)
            if x:
                ax.set_xlabel(x)
            if y:
                ax.set_ylabel(y)
            if xlim:
                ax.set_xlim(xlim)
            if ylim:
                ax.set_ylim(ylim)

            if keys:
                plt.key()

            if show:
                plt.show()
            return ax
        except ImportError:
            raise ImportError("require matplotlib")

    def append(self, dset, copy=True):
        """append anythin having `key` attribute
        """
        if copy:
            d0 = dset.copy()
        else:
            d0 = dset

        def check_key(self, key):
            for k0 in self.keys():
                if k0 == key:
                    raise KeyError("must have different key")

        if hasattr(dset, 'key'):
            check_key(self, dset.key)
            super(DatasetList, self).append(d0)
        elif isinstance(dset, dict):
            for key, values in iteritems(dset):
                check_key(self, key)
                super(DatasetList, self).append(DataArray({key: values}))
        else:
            raise ValueError()

    def remove(self, dset):
        for idx, d in enumerate(self):
            if dset.key == d.key:
                # do not work with
                # super(DatasetList, self).remove(d)
                # TypeError: 'NotImplementedType' object is not callable
                # why?
                super(DatasetList, self).remove(self.__getitem__(idx))

    def from_datasetlist(self, dslist, copy=True):
        self.from_sequence(dslist, copy=copy)
        return self

    def from_sequence(self, dslist, copy=True):
        for d in dslist:
            self.append(d, copy=copy)
        return self

    def chunk_average(self, n_chunks):
        return dict(map(lambda x: (x.key, x.chunk_average(n_chunks)), self))

    def topk(self, k):
        return dict((x.key, x.topk(k)) for x in self)

    def lowk(self, k):
        from heapq import nsmallest
        return dict((x.key, list(nsmallest(k, x))) for x in self)

    def head(self, k):
        return dict((x.key, x.head(k, restype='list')) for x in self)

    def tail(self, k):
        return dict((x.key, x.tail(k)) for x in self)

    def groupby(self, func_or_key):
        return _groupby(self, func_or_key)

    def cpptraj_dtypes(self):
        return np.array([x.cpptraj_dtype for x in self])

    def names(self):
        return np.array([x.name for x in self])

    def sort(self, key=None, copy=False, *args, **kwd):
        return from_sequence(sorted(self, key=key, *args, **kwd), copy=copy)

    def __add__(self, other):
        raise NotImplementedError()

    def __mul__(self, other):
        raise NotImplementedError()

    def extend(self, other, copy=True):
        self.from_sequence(other, copy=copy)

    def dot(self, idx0, idx1):
        """equal to np.dot(D[idx0].values, D[idx1].values)
        """
        return np.dot(self[idx0].values, self[idx1].values)

    def astype(self, dtype):
        return self.values.astype(dtype)

    def flatten(self):
        """return flatten numpy array
        """
        try:
            return self.values.flatten()
        except ValueError:
            from pytraj.tools import flatten
            import numpy as np
            return np.asarray(flatten(self))
