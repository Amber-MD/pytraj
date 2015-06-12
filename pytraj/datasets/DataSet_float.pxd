# distutils: language = c++
from libcpp.vector cimport vector
from .DataSet cimport _DataSet, DataSet
from .DataSet_1D cimport _DataSet_1D, DataSet_1D


cdef extern from "DataSet_float.h": 
    cdef cppclass _DataSet_float "DataSet_float" (_DataSet_1D):
        _DataSet_float() 
        @staticmethod
        _DataSet * Alloc() 
        float& operator[](size_t idx)
        float& index_opr "operator[]"(size_t idx)
        int Size()
        void Resize(size_t)

cdef class DataSet_float (DataSet_1D):
    cdef _DataSet_float* thisptr
    cdef bint py_free_mem 

