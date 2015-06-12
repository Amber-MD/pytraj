# distutils: language = c++
from libcpp.vector cimport vector
from .DataSet cimport _DataSet, DataSet
from .DataSet_1D cimport _DataSet_1D, DataSet_1D


cdef extern from "DataSet_double.h": 
    cdef cppclass _DataSet_double "DataSet_double" (_DataSet_1D):
        _DataSet_double() 
        @staticmethod
        _DataSet * Alloc() 
        double& operator[](size_t idx)
        double& index_opr "operator[]"(size_t idx)
        const vector[double]& Data() const 
        void assign_opr "operator =" (const vector[double]& rhs)
        void AddElement(double d)
        void Resize(size_t sizeIn)
        size_t Size()
        int Sync() 
        void Info() const 
        int Allocate1D(size_t)
        void Add(size_t, const void *)
        double Dval(size_t idx) const 
        double Xcrd(size_t idx) const 
        void Append(const _DataSet_double&)
        void SetNOE(double b, double bh, double r)
        double NOE_bound() const 
        double NOE_boundH() const 
        double NOE_rexp() const 
        void ShiftTorsions(double, double)


cdef class DataSet_double (DataSet_1D):
    cdef _DataSet_double* thisptr
    cdef bint py_free_mem 
