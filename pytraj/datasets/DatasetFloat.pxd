# distutils: language = c++
from libcpp.vector cimport vector
from .DataSet cimport _DataSet, DataSet
from .DataSet_1D cimport _DataSet_1D, DataSet_1D


cdef extern from "DataSet_float.h": 
    cdef cppclass _DatasetFloat "DataSet_float" (_DataSet_1D):
        _DatasetFloat() 
        @staticmethod
        _DataSet * Alloc() 
        float& operator[](size_t idx)
        float& index_opr "operator[]"(size_t idx)
        int Size()
        void Resize(size_t)

cdef class DatasetFloat (DataSet_1D):
    cdef _DatasetFloat* thisptr
    cdef bint py_free_mem 

