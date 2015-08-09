# distutils: language = c++
from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as incr
from cpython.array cimport array

cdef class DatasetMatrix3x3 (DataSet_1D):
    def __cinit__(self):
        # TODO : Use only one pointer? 
        self.baseptr0 = <_DataSet*> new _DatasetMatrix3x3()
        # make sure 3 pointers pointing to the same address?
        self.baseptr_1 = <_DataSet_1D*> self.baseptr0
        self.thisptr = <_DatasetMatrix3x3*> self.baseptr0

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
        if self.size <= 0:
            raise ValueError("size should be > 0")

        cdef Matrix_3x3 mat = Matrix_3x3()
        mat.thisptr[0] = self.thisptr[0][idx]
        return mat

    def __setitem__(self, int idx, double value):
        raise NotImplementedError()
        
    def __iter__(self):
        """return copy"""
        if self.size <= 0:
            raise ValueError("size should be > 0")
        cdef vector[_Matrix_3x3].iterator it = self.thisptr.begin()
        cdef Matrix_3x3 mat

        while it != self.thisptr.end(): 
            mat = Matrix_3x3()
            mat.thisptr[0] = deref(it)
            incr(it)
            yield mat

    def append(self, Matrix_3x3 mat):
        self.thisptr.AddMat3x3(mat.thisptr[0])

    def tolist(self):
        return self.to_ndarray().tolist()

    def to_pyarray(self):
        """slow"""
        return array('d', self.to_ndarray().flatten())

    def to_ndarray(self, copy=True):
        """return a copy
        """
        import numpy as np
        try:
            return np.array([x.to_ndarray(copy=copy) for x in self])
        except ValueError:
            return np.array([], dtype='f8')
        
