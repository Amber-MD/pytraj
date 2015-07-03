# distutils: language = c++
from __future__ import absolute_import
from .DataSet cimport _DataSet, DataSet, DataType


cdef extern from "DataSet_2D.h": 
    # DataSet_2D.h
    ctypedef enum MatrixType "DataSet_2D::MatrixType":
        pass
    ctypedef enum MatrixKind "DataSet_2D::MatrixKind":
        pass
    cdef cppclass _DataSet_2D "DataSet_2D" (_DataSet):
        _DataSet_2D() 
        _DataSet_2D(DataType tIn, int wIn, int pIn)
        # virtual methods
        int Allocate2D(size_t, size_t) 
        int AllocateHalf(size_t) 
        int AllocateTriangle(size_t) 
        double GetElement(size_t, size_t) const  
        size_t Nrows() const  
        size_t Ncols() const  
        double * MatrixArray() const  
        MatrixKind Kind "MatrixKind"() const  
        # end virtual methods

        void Add(size_t, const void *)
        const char * MatrixTypeString(MatrixType m)
        const char * MatrixOutputString(MatrixType m)

cdef class DataSet_2D (DataSet):
    cdef _DataSet_2D* baseptr_1
