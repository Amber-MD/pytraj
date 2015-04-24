# distutils: language = c++
from libcpp.vector cimport vector
from .DataSet cimport _DataSet, DataSet
from .DataSet_1D cimport _DataSet_1D, DataSet_1D
from ..CpptrajFile cimport _CpptrajFile, CpptrajFile


cdef extern from "DataSet_integer.h": 
    cdef cppclass _DataSet_integer "DataSet_integer" (_DataSet_1D):
        _DataSet_integer() 
        @staticmethod
        _DataSet * Alloc() 
        int& operator[](size_t idx)
        int& index_opr "operator[]"(size_t idx)
        void AddElement(int i)
        int Size()
        void Add( size_t, const void* )

cdef class DataSet_integer (DataSet_1D):
    cdef _DataSet_integer* thisptr
    cdef bint py_free_mem 

