# distutils: language = c++
from cpython.array cimport array as pyarray
from cython.view cimport array as cyarray

# python level
#from pytraj.optional_libs import HAS_NUMPY, ndarray

cdef class DataSet_float (DataSet_1D):
    def __cinit__(self):
        # TODO : Use only one pointer? 
        self.baseptr0 = <_DataSet*> new _DataSet_float()
        # make sure 3 pointers pointing to the same address?
        self.baseptr_1 = <_DataSet_1D*> self.baseptr0
        self.thisptr = <_DataSet_float*> self.baseptr0

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
        #return self.thisptr.index_opr(idx)
        return self.data[idx]

    def __setitem__(self, idx, value):
        self.data[idx] = value
        
    def __iter__(self):
        cdef int i
        for i in range(self.size):
            yield self.thisptr.index_opr(i)

    def resize(self, size_t sizeIn):
        self.thisptr.Resize(sizeIn)

    property data:
        def __get__(self):
            """return memoryview of data array
            """
            cdef cyarray myview
            cdef int size = self.size
            cdef float* ptr

            if size == 0:
                return None
            ptr = &self.thisptr.index_opr(0)
            myview = <float[:size]> ptr
            return myview

        def __set__(self, data):
            raise NotImplementedError()

    def append(self, ds):
        cdef int new_size = self.size + ds.size
        cdef int j
        self.resize(new_size)

        j = 0
        for i in range(self.size, new_size):
            self[i] = ds[j]
            j += 1
