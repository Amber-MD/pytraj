# distutils: language = c++
from libcpp.vector cimport vector
from .DataSet_2D cimport *
from ..math.Matrix cimport *

#ctypedef Matrix[double].iterator iterator
ctypedef vector[double] Darray

cdef extern from "DataSet_MatrixFlt.h": 
    cdef cppclass _DatasetMatrixFloat  "DataSet_MatrixFlt" (_DataSet_2D):
        _DataSet_MatrixFlt() 
        @staticmethod
        _DataSet * Alloc() 


cdef class DatasetMatrixFloat(DataSet_2D):
    cdef _DatasetMatrixFloat * thisptr
    cdef bint py_free_mem
