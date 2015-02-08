# distutils: language = c++
from libcpp.vector cimport vector
from pytraj.datasets.DataSet cimport _DataSet, DataSet
from pytraj.datasets.DataSet_1D cimport _DataSet_1D, DataSet_1D
from pytraj.CpptrajFile cimport _CpptrajFile, CpptrajFile


cdef extern from "DataSet_float.h": 
    cdef cppclass _DataSet_float "DataSet_float" (_DataSet_1D):
        _DataSet_float() 
        @staticmethod
        _DataSet * Alloc() 
        float& operator[](size_t idx)
        float& index_opr "operator[]"(size_t idx)
        int Size()


cdef class DataSet_float (DataSet_1D):
    cdef _DataSet_float* thisptr
    cdef bint py_free_mem 

