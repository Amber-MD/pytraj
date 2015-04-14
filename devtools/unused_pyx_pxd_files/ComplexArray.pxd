# distutils: language = c++
from pytraj.ArrayIterator cimport *


cdef extern from "ComplexArray.h": 
    cdef cppclass _ComplexArray "ComplexArray":
        _ComplexArray() 
        _ComplexArray(int)
        _ComplexArray(const _ComplexArray&)
        #_ComplexArray& operator =(const _ComplexArray&)
        #~_ComplexArray() 
        void Allocate(int)
        void Assign(const _ComplexArray&)
        void PadWithZero(int)
        void Normalize(double)
        void SquareModulus() 
        void ComplexConjTimes(const _ComplexArray&)
        double * CAptr() 
        int size() const 
        double& operator[](int idx)
        const double& operator[](int idx) const 
        #const iterator begin() const 
        #const iterator end() const 

cdef class ComplexArray:
    cdef _ComplexArray* thisptr
