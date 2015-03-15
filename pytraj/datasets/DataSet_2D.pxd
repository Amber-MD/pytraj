# distutils: language = c++
from pytraj.datasets.DataSet cimport _DataSet, DataSet, DataType
from pytraj.CpptrajFile cimport _CpptrajFile, CpptrajFile


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
        void Write2D(_CpptrajFile&, int, int) const  
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
