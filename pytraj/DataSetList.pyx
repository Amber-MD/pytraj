# distutils: language = c++
from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as incr
from cpython.array cimport array

# python level
from pytraj.datasets.cast_dataset import cast_dataset
from pytraj.utils.check_and_assert import _import, is_array
from pytraj.utils.check_and_assert import _import_numpy
from pytraj.utils.check_and_assert import is_word_in_class_name
from collections import defaultdict
from pytraj._utils cimport get_positive_idx
from pytraj.externals.six import string_types
from pytraj.compat import set
from pytraj.utils import is_int
from pytraj.exceptions import *
from pytraj.DataFile import DataFile
from pytraj.ArgList import ArgList

# can not import cpptraj_dict here
# if doing this, we introduce circle-import since cpptraj_dict already imported
# DataSet
from pytraj.cpptraj_dict import DataTypeDict

_, np = _import_numpy()

cdef class DataSetList:
    def __cinit__(self, py_free_mem=True):
        # py_free_mem is a flag to tell pytraj should free memory or let 
        # cpptraj does
        # check ./CpptrajState.pyx
        self.thisptr = new _DataSetList()
        self.py_free_mem = py_free_mem

    def __dealloc__(self):
        if self.py_free_mem:
            del self.thisptr

    def __str__(self):
        has_pd, _ = _import("pandas")
        if not has_pd:
            msg = "<pytraj.DataSetList with %s datasets> (install pandas for pretty print)" % self.size
            return msg
        else:
            df = self.to_dataframe()
            return df.__str__()

    def __repr__(self):
        return self.__str__()

    def copy(self):
        cdef DataSetList dnew = DataSetList()
        for d in self:
            dnew._add_copy_of_set(d)
        return dnew

    def __call__(self, *args, **kwd):
        return self.groupby(*args, **kwd)

    def clear(self):
        self.thisptr.Clear()

    def __iadd__(self, DataSetList other):
        self.thisptr.addequal(other.thisptr[0])
        return self

    def __iter__(self):
        cdef const_iterator it
        cdef DataSet dset
        it = self.thisptr.begin()

        while it != self.thisptr.end():
            dset = DataSet()
            dset.baseptr0 = deref(it)
            yield cast_dataset(dset, dtype=dset.dtype)
            incr(it)

    def __len__(self):
        cdef const_iterator it
        cdef DataSet dset
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

    def remove_set(self, DataSet dset):
        self.thisptr.RemoveSet(dset.baseptr0)

    def __getitem__(self, idx):
        """return a DataSet instance
        Memory view is applied (which mean this new insance is just alias of self[idx])
        Should we use a copy instead?
        """
        cdef DataSet dset = DataSet()
        cdef DataSetList new_dslist
        cdef int start, stop, step
        cdef object _idx # _idx can be 'int' or 'string'

        if self.size == 0:
            raise ValueError("size = 0: can not index")

        if is_int(idx):
            _idx = get_positive_idx(idx, self.size)
            # get memoryview
            dset.baseptr0 = self.thisptr.index_opr(_idx)
            return cast_dataset(dset, dtype=dset.dtype)
        elif isinstance(idx, string_types):
             # return a list of datasets having idx as legend
             for d0 in self:
                 if d0.legend.upper() == idx.upper():
                     return d0
        elif isinstance(idx, slice):
            # return new view of `self`
            start, stop, step = idx.indices(self.size)
            new_dslist = DataSetList()
            new_dslist.py_free_mem = False # view
            for _idx in range(start, stop, step):
                new_dslist.add_existing_set(self[_idx])
            return new_dslist
        elif is_array(idx) or isinstance(idx, list):
            new_dslist = DataSetList()
            new_dslist.py_free_mem = False # view
            for _idx in idx: 
                new_dslist.add_existing_set(self[_idx])
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

    def set_ensemble_num(self,int i):
        self.thisptr.SetEnsembleNum(i)

    def allocate_sets(self,long int i):
        self.thisptr.AllocateSets(i)

    def set_precision_of_data_sets(self, string nameIn, int widthIn, int precisionIn):
        self.thisptr.SetPrecisionOfDataSets(nameIn, widthIn, precisionIn)

    def get_reference_frame(self, name):
        cdef DataSet dset = DataSet()
        # return a view of dataset
        dset.baseptr0 = self.thisptr.GetReferenceFrame(name.encode())
        return dset

    def get_dataset(self, idx=None, name=None, dtype=None):
        """
        return DataSet instance
        Input:
        =====
        name :: str, optional
        idx :: integer, optional
        """
        cdef DataSet dset = DataSet()

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
        cdef DataSetList dlist = DataSetList()
        dlist.thisptr[0] = self.thisptr.GetMultipleSets(s)
        return dlist

    def add_set(self, dtype=None, name="", default_name=""):
        # TODO: check cpptraj for this method
        cdef DataSet dset = DataSet()
        if dtype is None:
            raise ValueError("dtype must not be None")
        dtype = dtype.upper()
        name = name.encode()
        default_name = default_name.encode()
        dset.baseptr0 = self.thisptr.AddSet(DataTypeDict[dtype], name, default_name)
        return dset

    def add_existing_set(self, DataSet ds):
        self.thisptr.AddSet(ds.baseptr0)
        
    def add_setidx(self, DataType inType, string nameIn, int idx):
        cdef DataSet dset = DataSet()
        dset.baseptr0 = self.thisptr.AddSetIdx(inType, nameIn, idx)
        if not dset.baseptr0:
            raise MemoryError("Can not initialize pointer")
        return dset

    def _add_copy_of_set(self, DataSet dset):
        self.thisptr.AddCopyOfSet(dset.baseptr0)

    def add_set_aspect(self, dtype, name=None, aspect=None):
        """add new dataset
        Paramters
        --------
        dtype : str
            DataType

        name_1 : str

        name_2 : str
        """
        cdef DataSet ds = DataSet()
        if aspect is None:
            aspect = name
        ds.baseptr0 = self.thisptr.AddSetAspect(DataTypeDict[dtype], 
                                                name.encode(), aspect.encode())
        return ds

    def find_coords_set(self, name):
        name = name.encode()
        cdef DataSet dset = DataSet()
        dset.baseptr0 = self.thisptr.FindCoordsSet(name)
        return dset

    def find_set_of_type(self, filename, dtype):
        cdef DataSet dset = DataSet()

        dtype = dtype.upper()
        dset.baseptr0 = self.thisptr.FindSetOfType(filename.encode(), DataTypeDict[dtype])

        if not dset.baseptr0:
            raise MemoryError("Can not initialize pointer")
        return dset

    # TODO: combine those methods into one
    def get_legends(self):
        """return a list"""
        tmp_list = []
        for d0 in self:
            tmp_list.append(d0.legend)
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

    def keys(self):
        return self.get_legends()

    def groupby(self, key, mode='legend'):
        """"return a new DataSetList object as a view of `self`

        Parameters
        ----------
        key : str
            keyword for searching
        mode: str, default='legend'
            mode = 'legend' | 'name' | 'dtype' | 'aspect'
        """
        import re
        cdef DataSetList dtmp

        dtmp = DataSetList()

        # dont free mem here
        dtmp.py_free_mem = False
        for d0 in self:
            att = getattr(d0, mode)
            if re.search(key, att):
                dtmp.add_existing_set(d0)
        return dtmp

    def tolist(self):
        """return a list of list/array"""
        try:
            return [d0.tolist() for d0 in self]
        except:
            raise PytrajConvertError("dont know how to convert to list")

    def to_dict(self, use_numpy=False):
        """return a dict object with key=legend, value=list"""
        try:
            if use_numpy:
                return dict((d0.legend, d0.to_ndarray()) for d0 in self)
            else:
                return dict((d0.legend, d0.tolist()) for d0 in self)
        except:
            raise PytrajConvertError("don't know tho to convert to dict")

    def to_ndarray(self):
        has_np, np = _import("numpy")
        if has_np:
            try:
                if self.size == 1:
                    return self[0].to_ndarray()
                else:
                    # more than one set
                    return np.asarray([d0.to_ndarray() for d0 in self])
            except:
                raise PytrajConvertError("don't know how to convert to ndarray")
        else:
            raise PytrajConvertError("don't have numpy")

    def to_dataframe(self):
        """return pandas' DataFrame

        Requires
        --------
        pandas
        """
        _, pandas = _import("pandas")
        my_dict = dict((d0.legend, d0.to_ndarray()) for d0 in self)
        return pandas.DataFrame(my_dict)

    def set_py_free_mem(self, bint value):
        # we only expose py_free_mem in cython (not pure python)
        # we don't want to change *.pxd signature files since this 
        # requires recompiling *pyx codes
        self.py_free_mem = value

    def _base_dataset_iter(self):
        """return a list of baseclass DataSet"""
        cdef const_iterator it
        cdef DataSet dset
        it = self.thisptr.begin()


        while it != self.thisptr.end():
            dset = DataSet()
            dset.baseptr0 = deref(it)
            yield dset
            incr(it)

    def apply(self, func):
        for d in self:
            arr = np.asarray(d.data)
            arr[:] = func(arr)

    def data(self):
        """"""
        raise NotImplementedError("Not yet")

    def mean(self, axis=1):
        return self.to_ndarray().mean(axis=axis)

    def median(self, axis=1):
        return np.median(self.to_ndarray(), axis=axis)

    def std(self, axis=1):
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

    def sum(self, legend=None, axis=1):
        _, np = _import_numpy()
        if not legend:
            return np.sum(self.to_ndarray(), axis=axis)
        else:
            return self.groupby(legend).sum(axis=axis)

    def cumsum(self, axis=1):
        """Return the cumulative sum of the elements along a given axis.
        (from numpy doc)
        """
        return np.cumsum(self.to_ndarray(), axis=axis)

    def count(self, number=None):
        return dict((d0.legend, d0.count(number)) for d0 in self)

    def read_data(self, filename, arg=""):
        df = DataFile()
        df.read_data(filename, ArgList(arg), self)
