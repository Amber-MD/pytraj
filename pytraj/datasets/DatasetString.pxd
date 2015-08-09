# distutils: language = c++
from libcpp.string cimport string
from .DataSet cimport _DataSet, DataSet
from .DataSet_1D cimport _DataSet_1D, DataSet_1D


cdef extern from "DataSet_string.h": 
    cdef cppclass _DatasetString "DataSet_string" (_DataSet_1D):
        _DatasetString()
        _DataSet * Alloc() 
        string& index_opr "operator[]"(size_t idx)
        void AddElement(const string& s)
        void Resize(size_t sizeIn)
        int Size()

cdef class DatasetString(DataSet_1D):
    cdef _DatasetString* thisptr
    cdef bint py_free_mem
