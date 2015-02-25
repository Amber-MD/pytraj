# distutils: language = c++
from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as incr

# python level
from pytraj.cast_dataset import cast_dataset
from pytraj.utils.check_and_assert import _import
from collections import defaultdict

# can not import cpptraj_dict here
# if doing this, we introduce circle-import since cpptraj_dict already imported
# DataSet
from pytraj.cpptraj_dict import DataTypeDict

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

    def __cal__(self, *args, **kwd):
        return self.get_dataset(*args, **kwd)

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

    property size:
        def __get__(self):
            return self.thisptr.size()

    def ensemble_num(self):
        return self.thisptr.EnsembleNum()

    def remove_set(self, DataSet dset):
        self.thisptr.RemoveSet(dset.baseptr0)

    def __getitem__(self, int idx):
        """return a DataSet instance
        Memory view is applied (which mean this new insance is just alias of self[idx])
        Should we use a copy instead?
        """
        cdef DataSet dset = DataSet()
        if idx >= len(self) or idx < 0:
            raise ValueError("index is out of range")
        # get memoryview
        dset.baseptr0 = self.thisptr.index_opr(idx)
        return cast_dataset(dset, dtype=dset.dtype)

    def set_ensemble_num(self,int i):
        self.thisptr.SetEnsembleNum(i)

    def allocate_sets(self,long int i):
        self.thisptr.AllocateSets(i)

    def set_precision_of_data_sets(self, string nameIn, int widthIn, int precisionIn):
        self.thisptr.SetPrecisionOfDataSets(nameIn, widthIn, precisionIn)

    def get_set(self, string dsname, int idx, string attr_arg):
        cdef DataSet dset = DataSet()
        dset.baseptr0 = self.thisptr.GetSet(dsname, idx, attr_arg)

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
                    if d0.dtype == dtype:
                        dlist.append(d0[:])
                # return a list of arrays
                has_numpy, np = _import('numpy')
                if has_numpy:
                    return np.array(dlist)
                else:
                    return dlist

    def get_multiple_sets(self, string s):
        """TODO: double-check cpptraj"""
        cdef DataSetList dlist = DataSetList()
        dlist.thisptr[0] = self.thisptr.GetMultipleSets(s)
        return dlist

    def add_set(self, DataType dtype, string s, char * c):
        # TODO: check cpptraj for this method
        cdef DataSet dset = DataSet()
        dset.baseptr0 = self.thisptr.AddSet(dtype, s, c)
        return dset
        
    def add_setidx(self, DataType inType, string nameIn, int idx):
        cdef DataSet dset = DataSet()
        dset.baseptr0 = self.thisptr.AddSetIdx(inType, nameIn, idx)
        if not dset.baseptr0:
            raise MemoryError("Can not initialize pointer")
        return dset

    def add_copy_of_set(self, DataSet dset):
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
        ds.baseptr0 = self.thisptr.AddSetAspect(DataTypeDict[dtype], name.encode(), aspect.encode())
        return ds

    def printlist(self):
        self.thisptr.List()

    def find_coords_set(self, string filename):
        cdef DataSet dset = DataSet()
        dset.baseptr0 = self.thisptr.FindCoordsSet(filename)
        if not dset.baseptr0:
            raise MemoryError("Can not initialize pointer")
        return dset
