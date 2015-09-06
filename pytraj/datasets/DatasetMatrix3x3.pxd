# distutils: language = c++
from libcpp.vector cimport vector
from .base cimport _DataSet, DataSet
from .dataset_1d cimport _DataSet_1D, DataSet_1D
from ..math.Matrix_3x3 cimport _Matrix_3x3, Matrix_3x3


cdef extern from "DataSet_Mat3x3.h": 
    ctypedef vector[_Matrix_3x3].iterator mat_iterator
    cdef cppclass _DatasetMatrix3x3 "DataSet_Mat3x3" (_DataSet_1D):
        _DatasetMatrix3x3()
        @staticmethod
        _DataSet * Alloc() 
        bint Empty()
        void AddMat3x3(_Matrix_3x3)
        mat_iterator begin()
        mat_iterator end()
        _Matrix_3x3& operator[](int i)


cdef class DatasetMatrix3x3(DataSet_1D):
    cdef _DatasetMatrix3x3* thisptr
    cdef bint py_free_mem 
