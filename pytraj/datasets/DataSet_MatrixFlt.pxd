# distutils: language = c++
from libcpp.vector cimport vector
from pytraj.datasets.DataSet_2D cimport *
from pytraj.Matrix cimport *

#ctypedef Matrix[double].iterator iterator
ctypedef vector[double] Darray

cdef extern from "DataSet_MatrixFlt.h": 
    cdef cppclass _DataSet_MatrixFlt "DataSet_MatrixFlt" (_DataSet_2D):
        _DataSet_MatrixFlt() 
        @staticmethod
        _DataSet * Alloc() 


cdef class DataSet_MatrixFlt (DataSet_2D):
    cdef _DataSet_MatrixFlt* thisptr
    cdef bint py_free_mem
