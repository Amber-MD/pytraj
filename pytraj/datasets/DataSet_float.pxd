# distutils: language = c++
from DataSet_1D cimport *


cdef extern from "DataSet_float.h": 
    cdef cppclass _DataSet_float "DataSet_float":
        _DataSet_float() : _DataSet_1D(FLOAT, 8, 3)
        _DataSet * Alloc() 
        float & operator [ ](size_t idx)
        float operator [ ](size_t idx)const 
        void AddElement(float f)
        void Resize(size_t sizeIn)
        size_t Size() const 
        int Sync() 
        void Info() const 
        int Allocate1D(size_t)
        void Add(size_t, const void *)
        double Dval(size_t idx)const 
        double Xcrd(size_t idx)const 
        void WriteBuffer(_CpptrajFile &, size_t)const 
