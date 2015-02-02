# distutils: language = c++
from DataSet_2D cimport *
from Matrix cimport *


cdef extern from "DataSet_MatrixFlt.h": 
    cdef cppclass _DataSet_MatrixFlt "DataSet_MatrixFlt":
        _DataSet_MatrixFlt() : _DataSet_2D(MATRIX_FLT, 12, 4)
        float & operator [ ](size_t idx)
        _DataSet * Alloc() 
        size_t Size() const 
        int Sync() 
        void Info() const 
        int Allocate2D(size_t x, size_t y)
        int AllocateHalf(size_t x)
        int AllocateTriangle(size_t x)
        void Write2D(_CpptrajFile &, int, int)const 
        double GetElement(size_t x, size_t y)const 
        size_t Nrows() const 
        size_t Ncols() const 
        double * MatrixArray() const 
        MatrixKind Kind() const 
        MatrixType Type() const 
        int AddElement(float d)
        void SetElement(size_t x, size_t y, float d)
