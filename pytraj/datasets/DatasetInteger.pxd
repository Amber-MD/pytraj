# distutils: language = c++
from libcpp.vector cimport vector
from .DataSet cimport _DataSet, DataSet
from .DataSet_1D cimport _DataSet_1D, DataSet_1D


cdef extern from "DataSet_integer.h": 
    cdef cppclass _DatasetInteger "DataSet_integer" (_DataSet_1D):
        _DatasetInteger() 
        @staticmethod
        _DataSet * Alloc() 
        int& operator[](size_t idx)
        int& index_opr "operator[]"(size_t idx)
        void AddElement(int i)
        int Size()
        void Resize(size_t)
        void Add( size_t, const void* )

cdef class DatasetInteger (DataSet_1D):
    cdef _DatasetInteger* thisptr
    cdef bint py_free_mem 

