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

    def __getitem__(self, int idx):
        #return self.thisptr.index_opr(idx)
        # use self.data so we can use fancy indexing
        return self.data[idx].decode()

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
    def size(self):
        return self.thisptr.Size()

