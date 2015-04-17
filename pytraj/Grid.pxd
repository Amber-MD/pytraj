# distutil: language = c++

from pytraj.ArrayIterator cimport *

cdef extern from "Grid.h":
    #ctypedef GridIterator "Grid::ArrayIterator[T]" iterator
    cdef cppclass _Grid "Grid" [T]:
        _Grid()
        #~Grid() 
        _Grid(const _Grid &)
        #not supported yet>
        #Grid & operator =(const Grid &)
        #T & operator [](size_t idx)
        T& index_opr "operator[]"(size_t idx)
        T& operator[](size_t idx)
        size_t size() const 
        int resize(size_t, size_t, size_t)
        size_t NX() const 
        size_t NY() const 
        size_t NZ() const 
        long int incrementBy(int, int, int, const T &)
        void setGrid(int, int, int, const T &)
        const T& element(int, int, int)const 
        long int CalcIndex(int x, int y, int z)const 
        #iterator begin() 
        #iterator end() 

cdef class Grid:
    pass
    cdef _Grid[float]* thisptr
