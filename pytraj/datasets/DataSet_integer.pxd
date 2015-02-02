# distutils: language = c++
from DataSet_1D cimport *


cdef extern from "DataSet_integer.h": 
    cdef cppclass _DataSet_integer "DataSet_integer":
        _DataSet_integer()
        _DataSet * Alloc() 
        int & operator [ ](size_t idx)
        int operator [ ](size_t idx)const 
        void AddElement(int i)
        void Resize(size_t sizeIn)
        inline void AddVal(size_t, int)
        size_t Size() const 
        int Sync() 
        void Info() const 
        int Allocate1D(size_t)
        void Add(size_t, const void *)
        double Dval(size_t idx)const 
        double Xcrd(size_t idx)const 
        void WriteBuffer(_CpptrajFile &, size_t)const 
        #iterator begin() 
        #iterator end() 
