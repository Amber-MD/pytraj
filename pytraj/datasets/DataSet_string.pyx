# distutils: language = c++
from cpython.array cimport array as pyarray

# python level
#from pytraj.optional_libs import HAS_NUMPY, ndarray

cdef class DataSet_string (DataSet_1D):
    def __cinit__(self):
        self.baseptr0 = <_DataSet*> new _DataSet_string()
        # make sure 3 pointers pointing to the same address?
        self.baseptr_1 = <_DataSet_1D*> self.baseptr0
        self.thisptr = <_DataSet_string*> self.baseptr0

        # let Python/Cython free memory
        self.py_free_mem = True

    def __dealloc__(self):
        if self.py_free_mem:
            del self.thisptr

    def alloc(self):
        '''return a memoryview as DataSet instane'''
        cdef DataSet dset = DataSet()
        dset.baseptr0 = self.thisptr.Alloc()
        return dset

    def __getitem__(self, idx):
        # this dataset does not have `data` attribute
        cdef int _idx
        cdef list slist = []

        if isinstance(idx, slice):
            if (idx.start, idx.stop, idx.step) is (None, None, None):
                for _idx in self.size:
                    slist.append(self[_idx])
                return slist 
            else:
                raise IndexError("only support `:` for slice")
        elif isinstance(idx, (long, int)):
            return self.thisptr.index_opr(idx)
        else:
            raise IndexError("index must be integer or `:` ")

    def __setitem__(self, int idx, value):
        cdef string* ptr
        ptr = &(self.thisptr.index_opr(idx))
        ptr[0] = value
        
    def __iter__(self):
        cdef int i
        for i in range(self.size):
            yield self.thisptr.index_opr(i)

    def resize(self, size_t sizeIn):
        self.thisptr.Resize(sizeIn)

    @property
    def data(self):
        """return 1D python array of `self`"""
        cdef pyarray arr = pyarray('i', [])
        cdef int i

        for i in range(self.baseptr0.Size()):
            arr.append(self[i])
        return arr
