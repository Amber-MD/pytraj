# distutils: language = c++
from __future__ import absolute_import

from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as incr
from cpython.array cimport array
from ..trajs.TrajectoryCpptraj cimport TrajectoryCpptraj

# python level
import numpy as np
from collections import defaultdict

from .._cyutils import get_positive_idx
from collections import defaultdict
from .cast_dataset import cast_dataset
from ..utils.check_and_assert import is_array
from ..externals.six import string_types
from ..compat import set
from ..utils import is_int
from ..datafiles.datafiles import DataFile
from ..core.cpptraj_core import ArgList

from pytraj.cpptraj_dict import DataTypeDict

__all__ = ['DatasetList']


cdef class DatasetList:
    """
    DatasetList holds data from cpptraj

    Notes
    -----
    Methods require pandas:
        * describe
        * to_dataframe
    """
    def __cinit__(self, py_free_mem=True):
        # py_free_mem is a flag to tell pytraj should free memory or let 
        # cpptraj does
        # check ./CpptrajState.pyx
        self.thisptr = new _DatasetList()
        self.py_free_mem = py_free_mem
        # Not all DatasetLists own their own data (if it's a MemoryView) for
        # instance, so this allows us to keep references to parent objects to
        # prevent them from getting GCed while their memory is still being used.
        self._parent_lists = []

    def __dealloc__(self):
        if self.py_free_mem:
            del self.thisptr

    def __str__(self):
        safe_msg = "<pytraj.DatasetList with %s datasets>" % self.size
        if self.size == 0:
            return safe_msg
        msg = "<pytraj.datasets.DatasetList - %s datasets>" % self.size
        return msg

    def __repr__(self):
        return self.__str__()

    def __contains__(self, Dataset other):
        cdef Dataset d0
        for d0 in self:
            if d0.baseptr0 == other.baseptr0:
                return True
        return False

    def copy(self):
        cdef DatasetList dnew = DatasetList()
        for d in self:
            dnew._add_copy_of_set(d)
        return dnew

    def __call__(self, *args, **kwd):
        return self.groupby(*args, **kwd)

    def clear(self):
        self.thisptr.Clear()

    def __iadd__(self, DatasetList other):
        self.thisptr.addequal(other.thisptr[0])
        return self

    def __iter__(self):
        cdef const_iterator it
        cdef Dataset dset
        it = self.thisptr.begin()

        while it != self.thisptr.end():
            dset = Dataset()
            dset.baseptr0 = deref(it)
            try:
                yield cast_dataset(dset, dtype=dset.dtype)
            except NotImplementedError:
                yield dset
            incr(it)

    def __len__(self):
        cdef const_iterator it
        cdef Dataset dset
        cdef int i
        it = self.thisptr.begin()

        i = 0
        while it != self.thisptr.end():
            i += 1
            incr(it)
        return i

    def len(self):
        return self.__len__()
            
    def is_empty(self):
        return self.thisptr.empty()

    @property
    def size(self):
        return self.thisptr.size()

    def ensemble_num(self):
        return self.thisptr.EnsembleNum()

    def remove_set(self, Dataset dset):
        self.thisptr.RemoveSet(dset.baseptr0)

    def __getitem__(self, idx):
        """return a Dataset instance
        Memory view is applied (which mean this new insance is just alias of self[idx])
        Should we use a copy instead?
        """
        cdef Dataset dset = Dataset()
        cdef int start, stop, step
        cdef object _idx # _idx can be 'int' or 'string'

        if self.size == 0:
            raise ValueError("size = 0: can not index")

        if is_int(idx):
            _idx = get_positive_idx(idx, self.size)
            # get memoryview
            dset.baseptr0 = self.thisptr.index_opr(_idx)
            dtmp = cast_dataset(dset, dtype=dset.dtype)
            dtmp._base = self
            return dtmp
        elif isinstance(idx, string_types):
             # return a list of datasets having idx as legend
             for d0 in self:
                 if d0.key.upper() == idx.upper():
                     d0._base = self
                     return d0
        elif isinstance(idx, slice):
            # return new view of `self`
            start, stop, step = idx.indices(self.size)
            new_dslist = self.__class__()
            new_dslist.set_py_free_mem(False)
            for _idx in range(start, stop, step):
                new_dslist.add_existing_set(self[_idx])
            new_dslist._parent_lists_append(self)
            return new_dslist
        elif is_array(idx) or isinstance(idx, list):
            new_dslist = self.__class__()
            new_dslist.set_py_free_mem(False)
            for _idx in idx: 
                new_dslist.add_existing_set(self[_idx])
            new_dslist._parent_lists_append(self)
            return new_dslist
        else:
            raise ValueError()

    def __setitem__(self, idx, value):
        """data copy

        dslist[0] += 1.
        """
        if hasattr(value, '_npdata'):
            self[idx]._npdata[:] = value._npdata[:]
        else:
            self[idx]._npdata[:] = value

    def allocate_sets(self,long int i):
        self.thisptr.AllocateSets(i)

    def get_dataset(self, idx=None, name=None, dtype=None):
        """
        return Dataset instance
        Input:
        =====
        name :: str, optional
        idx :: integer, optional
        """
        cdef Dataset dset = Dataset()

        if name is not None and idx is not None:
            raise ValueError("name and idx must not be set at the same time")
        else:
            if dtype is None: 
                if name is not None:
                    name = name.encode()
                    dset.baseptr0 = self.thisptr.GetDataSet(name)
                if idx is not None:
                    dset.baseptr0 = self.thisptr.index_opr(idx)
                return dset
            else:
                assert idx == None
                assert name == None
                dtype = dtype.upper()
                dlist = []
                for d0 in self:
                    if d0.dtype.upper() == dtype:
                        dlist.append(d0[:])
                # return a list of arrays
                return dlist

    def get_multiple_sets(self, string s):
        """TODO: double-check cpptraj"""
        cdef DatasetList dlist = DatasetList()
        dlist.thisptr[0] = self.thisptr.GetMultipleSets(s)
        return dlist

    def add_set(self, dtype=None, name="", default_name=""):
        """create new (empty) Dataset and add to `self`
        this is for internal use
        """
        cdef Dataset dset = Dataset()
        if dtype is None:
            raise ValueError("dtype must not be None")
        dtype = dtype.upper()
        name = name.encode()
        default_name = default_name.encode()
        dset.baseptr0 = self.thisptr.AddSet(DataTypeDict[dtype], name, default_name)
        return dset

    def add_existing_set(self, Dataset ds):
        self.thisptr.AddSet(ds.baseptr0)

    def _add_copy_of_set(self, Dataset dset):
        self.thisptr.AddCopyOfSet(dset.baseptr0)

    def _add_traj(self, TrajectoryCpptraj traj):
        self.thisptr.AddCopyOfSet(<_Dataset*> traj.thisptr)

    def add_copy_of_set(self, Dataset dset):
        self.thisptr.AddCopyOfSet(dset.baseptr0)

    def find_coords_set(self, name):
        name = name.encode()
        cdef Dataset dset = Dataset()
        dset.baseptr0 = self.thisptr.FindCoordsSet(name)
        return dset

    def find_set_of_type(self, filename, dtype):
        cdef Dataset dset = Dataset()

        dtype = dtype.upper()
        dset.baseptr0 = self.thisptr.FindSetOfType(filename.encode(), DataTypeDict[dtype])

        if not dset.baseptr0:
            raise MemoryError("Can not initialize pointer")
        return dset

    def get_legends(self):
        """return a list"""
        tmp_list = []
        for d0 in self._base_dataset_iter():
            tmp_list.append(d0.key)
        return tmp_list

    def get_aspects(self, is_set=True):
        """return a set of uniqure aspects if "is_set" = True
        else: return a full list
        """

        tmp_list = []
        for d0 in self._base_dataset_iter():
            tmp_list.append(d0.aspect)
        if is_set:
            return set(tmp_list)
        else:
            return tmp_list

    def get_scalar_types(self):
        """return a list"""
        tmp_list = []
        for d0 in self._base_dataset_iter():
            tmp_list.append(d0.scalar_type)
        return tmp_list

    def get_scalar_modes(self):
        """return a list"""
        tmp_list = []
        for d0 in self._base_dataset_iter():
            tmp_list.append(d0.scalar_mode)
        return tmp_list

    def get_dtypes(self):
        """return a list"""
        tmp_list = []
        for d0 in self._base_dataset_iter():
            tmp_list.append(d0.dtype)
        return tmp_list

    def keys(self):
        return self.get_legends()

    def iteritems(self):
        from pytraj.compat import zip
        for key in self.keys():
            yield key, self[key]

    def map(self, func):
        for d0 in self:
            yield func(d0)

    def filter(self, func, *args, **kwd):
        """return a new view of DatasetList of func return True"""
        dslist = self.__class__()

        if isinstance(func, (string_types, list, tuple)):
            return self.grep(func, *args, **kwd)
        elif callable(func):
            for d0 in self:
                if func(d0, *args, **kwd):
                    dslist.append(d0)
            return dslist
        else:
            raise NotImplementedError("func must be a string or callable")

    def grep(self, key, mode='key'):
        """"return a new DatasetList object as a view of `self`. This method is mostly used for testing.

        Parameters
        ----------
        key : str or list
            keyword for searching
        mode: str, default='key'
            mode = 'legend' | 'name' | 'dtype' | 'aspect'
        """
        import re

        # use __class__ so we can `groupby` return the same class
        # we subclass this Cython class to python level
        dtmp = self.__class__()

        # dont free mem here
        dtmp.set_py_free_mem(False)
        for d0 in self._base_dataset_iter():
            att = getattr(d0, mode)
            if isinstance(key, string_types):
                if re.search(key, att):
                    dtmp.add_existing_set(d0)
            elif isinstance(key, (list, tuple)):
                for _key in key:
                    if re.search(_key, att):
                        dtmp.add_existing_set(d0)
            else:
                raise ValueError("support string or list/tuple of strings")

        # dtmp is just a view, so keep track of parent to avoid GC
        dtmp._parent_lists_append(self)

        return dtmp

    def tolist(self):
        """return a list of list/array"""
        try:
            return [d0.tolist() for d0 in self]
        except:
            raise ValueError("dont know how to convert to list")

    def to_dict(self, use_numpy=False):
        """return a dict object with key=legend, value=list"""
        from collections import OrderedDict as dict
        try:
            if use_numpy:
                return dict((d0.key, d0.to_ndarray(copy=True)) for d0 in self)
            else:
                return dict((d0.key, d0.tolist()) for d0 in self)
        except:
            raise ValueError("don't know tho to convert to dict")

    @property
    def values(self):
        """return read-only ndarray"""
        # read-only
        try:
            return self.to_ndarray()
        except ValueError:
            raise ValueError("don't know how to cast to numpy array"
                             "try `tolist`, `to_dict`")

    def to_ndarray(self):
        """
        Notes: require numpy
        """
        # make sure to use copy=True to avoid memory error for memoryview
        try:
            if self.size == 1:
                return self[0].to_ndarray(copy=True)
            else:
                # more than one set
                return np.asarray([d0.to_ndarray(copy=True) for d0 in self])
        except:
            raise ValueError("don't know how to convert to ndarray")

    def to_dataframe(self):
        """return pandas' DataFrame

        Requires
        --------
        pandas
        """
        import pandas
        from collections import OrderedDict
        my_dict = OrderedDict((d0.key, d0.to_ndarray(copy=True)) for d0 in self)
        return pandas.DataFrame(my_dict)

    def set_py_free_mem(self, bint value):
        # we only expose py_free_mem in cython (not pure python)
        # we don't want to change *.pxd signature files since this 
        # requires recompiling *pyx codes
        self.py_free_mem = value

    def _base_dataset_iter(self):
        """return a list of baseclass Dataset"""
        cdef const_iterator it
        cdef Dataset dset
        it = self.thisptr.begin()


        while it != self.thisptr.end():
            dset = Dataset()
            dset.baseptr0 = deref(it)
            yield dset
            incr(it)

    def apply(self, func):
        for d in self:
            arr = np.asarray(d.data)
            arr[:] = func(arr)

    def mean(self, axis=1):
        """
        Notes: require numpy
        """
        return self.to_ndarray().mean(axis=axis)

    def median(self, axis=1):
        """
        Notes: require numpy
        """
        return np.median(self.to_ndarray(), axis=axis)

    def std(self, axis=1):
        """
        Notes: require numpy
        """
        return np.std(self.to_ndarray(), axis=axis)

    def min(self):
        arr = array('d', [])
        for d in self:
            arr.append(d.min())
        return arr

    def max(self):
        arr = array('d', [])
        for d in self:
            arr.append(d.max())
        return arr

    def sum(self, key=None, axis=1):
        """
        Notes: require numpy
        """
        if not key:
            return np.sum(self.to_ndarray(), axis=axis)
        else:
            return self.groupby(key).sum(axis=axis)

    def cumsum(self, axis=1):
        """Return the cumulative sum of the elements along a given axis.
        (from numpy doc)
        """
        return np.cumsum(self.to_ndarray(), axis=axis)

    def mean_with_error(self, other):

        ddict = defaultdict(tuple)
        for key, dset in self.iteritems():
            ddict[key] = dset.mean_with_error(other[key])
        return ddict

    def count(self, number=None):
        return dict((d0.key, d0.count(number)) for d0 in self)

    def read_data(self, filename, arg=""):
        df = DataFile()
        df.read_data(filename, ArgList(arg), self)

    # pandas related
    def describe(self):
        import pandas as pd
        return self.to_dataframe().describe()

    def write_all_datafiles(self, filenames=None):
        from pytraj.core.DataFileList import DataFileList
        df = DataFileList()

        for idx, d in enumerate(self):
            if filenames is None:
                # make default name
                d.key= d.key.replace(":", "_")
                fname = "pytraj_datafile_" + d.key + ".txt"
            else:
                fname = filenames[idx]
            df.add_dataset(fname, d)
        df.write_all_datafiles()

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
        np.savetxt(filename, values, format=formats, header=headers) 

    def _parent_lists_append(self, data):
        self._parent_lists.append(data)

    def _parent_lists_free(self):
        self._parent_lists_free = []

    def to_pyarray(self):
        if self.size > 1:
            return [d.to_pyarray() for d in self]
        else:
            return self[0].to_pyarray()
