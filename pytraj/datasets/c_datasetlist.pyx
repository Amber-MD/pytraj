# distutils: language = c++
from __future__ import absolute_import

from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as incr
from cpython.array cimport array
from ..c_traj.c_trajectory cimport TrajectoryCpptraj

# python level
import numpy as np
from collections import defaultdict

from ..cyutils import get_positive_idx
from collections import defaultdict
from .cast_dataset import cast_dataset
from ..utils.check_and_assert import is_array
from ..externals.six import string_types
from ..compat import set
from ..utils import is_int
from ..datafiles.datafiles import DataFile
from ..core.c_core import ArgList

from pytraj.c_dict import DataTypeDict

__all__ = ['DatasetList']


cdef class DatasetList:
    """This class behave like a simple Python list (you
    can index it by number) and a Python dictionary (you can index it by key). It's quite simple and
    is introduced here to hold data from cpptraj. Mostly for internal use.

    Examples
    --------
    >>> import pytraj as pt
    >>> from pytraj.datasets import CpptrajDatasetList
    >>> # create CpptrajDatasetList to store data
    >>> dslist = CpptrajDatasetList()
    >>> actlist = pt.ActionList(['distance :2 :3', 'vector :3 :5'], top=traj.top, dslist=dslist)
    >>> for frame in traj: actlist.do_actions(frame)
    >>> print(dslist)
    >>> print(dslist[0])
    """

    def __cinit__(self, _own_memory=True):
        # _own_memory is a flag to tell pytraj should free memory or let
        # cpptraj does
        # check ./CpptrajState.pyx
        self.thisptr = new _DatasetList()
        self._own_memory = _own_memory
        # Not all DatasetLists own their own data (if it's a MemoryView) for
        # instance, so this allows us to keep references to parent objects to
        # prevent them from getting GCed while their memory is still being used.
        self._parent_lists = []

    def __dealloc__(self):
        if self._own_memory:
            del self.thisptr

    def __str__(self):
        msg = "<pytraj.datasets.CpptrajDatasetList - %s datasets>" % self.size
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

    def clear(self):
        self.thisptr.Clear()

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

    def is_empty(self):
        return self.thisptr.empty()

    @property
    def size(self):
        return self.thisptr.size()

    def remove_set(self, Dataset dset):
        self.thisptr.RemoveSet(dset.baseptr0)

    def _pop(self, int i):
        self.remove_set(self[i])

    def __getitem__(self, idx):
        """return a Dataset instance
        Memory view is applied (which mean this new insance is just alias of self[idx])
        Should we use a copy instead?
        """
        cdef Dataset dset = Dataset()
        cdef int start, stop, step
        cdef object _idx  # _idx can be 'int' or 'string'

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
            new_dslist.set_own_memory(False)
            for _idx in range(start, stop, step):
                new_dslist.add_existing_set(self[_idx])
            new_dslist._parent_lists_append(self)
            return new_dslist
        elif is_array(idx) or isinstance(idx, list):
            new_dslist = self.__class__()
            new_dslist.set_own_memory(False)
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
                assert idx is None
                assert name is None
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

    def add_new(self, dtype=None, name=""):
        """create new (empty) Dataset and add to `self`

        Examples
        --------
        >>> from pytraj.datasets import CpptrajDatasetList
        >>> dslist = CpptrajDatasetList()
        >>> # add DatasetTopoloy
        >>> d0 = dslist.add_new('topology', name='myparm')
        >>> # add new Topology to d0
        >>> d0.data = traj.top
        >>> print(dslist[0])
        <pytraj.datasets.DatasetTopology: size=5293, key=myparm>
        """
        default_name = ""
        cdef Dataset dset = Dataset()
        if dtype is None:
            raise ValueError("dtype must not be None")
        dtype = dtype.upper()
        name = name.encode()
        default_name = default_name.encode()
        dset.baseptr0 = self.thisptr.AddSet(DataTypeDict[dtype], name, default_name)
        return cast_dataset(dset, dtype=dset.dtype)

    def add(self, *args, **kwd):
        '''alias of self.add_new
        '''
        return self.add_new(*args, **kwd)

    def add_existing_set(self, Dataset ds):
        self.thisptr.AddSet(ds.baseptr0)

    def _add_copy_of_set(self, Dataset dset):
        self.thisptr.AddCopyOfSet(dset.baseptr0)

    def _add_traj(self, TrajectoryCpptraj traj, name='_traj'):
        '''add an existing TrajectoryCpptraj with given name
        '''
        cdef _MetaData metadata
        metadata.SetName(name.encode())
        traj._own_memory = False
        (<_Dataset*> traj.thisptr).SetMeta(metadata)
        # self.thisptr.AddCopyOfSet(<_Dataset*> traj.thisptr)
        self.thisptr.AddSet(<_Dataset*> traj.thisptr)

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

    def filter(self, func, *args, **kwd):
        """return a new view of DatasetList of func return True"""
        dslist = self.__class__()
        dslist.set_own_memory(False)

        if isinstance(func, (string_types, list, tuple)):
            dslist = self.grep(func, *args, **kwd)
        elif callable(func):
            for d0 in self:
                if func(d0, *args, **kwd):
                    dslist.add_existing_set(d0)
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
        dtmp.set_own_memory(False)
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

    def to_dict(self):
        """return a dict object with key=legend, value=list"""
        from collections import OrderedDict
        try:
            return OrderedDict((d0.key, d0.to_ndarray(copy=True)) for d0 in self)
        except ValueError:
            try:
                # we use values (might be copy of dataset's iternal data)
                return OrderedDict((d0.key, d0.values) for d0 in self)
            except ValueError:
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

    def set_own_memory(self, bint value):
        # we only expose _own_memory in cython (not pure python)
        # we don't want to change *.pxd signature files since this
        # requires recompiling *pyx codes
        self._own_memory = value

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

    def read_data(self, filename, arg=""):
        df = DataFile()
        df.read_data(filename, ArgList(arg), self)

    def _parent_lists_append(self, data):
        self._parent_lists.append(data)

    def _parent_lists_free(self):
        self._parent_lists_free = []
