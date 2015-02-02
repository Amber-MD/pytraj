# distutil: language = c++
from pytraj.ArrayIterator cimport *

cdef extern from "Matrix.h":
    cdef cppclass _Matrix "Matrix" [T]:
        _Matrix()
        _Matrix(const _Matrix&)
        T& operator[](size_t idx)
        const T& operator[] (size_t)
        size_t size()
        int resize(size_t, size_t)
        const T& element(int, int) const
        T& element(int, int)
        size_t Nrows()
        size_t Ncols()
        int addElement(const T&)
        void setElement(int, int, const T&)
        const T* Ptr()
        T* Ptr()
        size_t CalcIndex(int, int)
        #Cython has not yet support nested ctypedef
        #ctypedef ArrayIterator[T] iterator
        #ArrayIterator[T] begin()
        #ArrayIterator[T] end()
