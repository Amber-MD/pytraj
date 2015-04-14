# distutils: language = c++
from libcpp.vector cimport vector
from pytraj.DataSetList cimport *
from pytraj.datasets.DataSet_1D cimport *
from pytraj.ArgList cimport *


cdef extern from "Array1D.h": 
    cdef cppclass _Array1D "Array1D":
        _Array1D() 
        _Array1D(const _Array1D&)
        #_Array1D& operator =(const _Array1D&)
        _Array1D(const _DataSetList&)
        size_t DetermineMax() const 
        int push_back(_DataSet_1D * const&)
        #_DataSet_1D * const& operator[](int idx) const 
        _DataSet_1D* index_opr "operator[]"(int idx)
        bint empty() const 
        #const_iterator begin() const 
        #const_iterator end() const 
        size_t size() const 
        void clear() 
        void SortArray1D() 
        int AddDataSets(const _DataSetList&)
        int AddTorsionSets(const _DataSetList&)
        int AddSetsFromArgs(const _ArgList&, const _DataSetList&)
        int CheckXDimension() const 


cdef class Array1D:
    cdef _Array1D* thisptr
