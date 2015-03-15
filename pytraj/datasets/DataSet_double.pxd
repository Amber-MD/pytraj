# distutils: language = c++
from libcpp.vector cimport vector
from pytraj.datasets.DataSet cimport _DataSet, DataSet
from pytraj.datasets.DataSet_1D cimport _DataSet_1D, DataSet_1D
from pytraj.CpptrajFile cimport _CpptrajFile, CpptrajFile


cdef extern from "DataSet_double.h": 
    cdef cppclass _DataSet_double "DataSet_double" (_DataSet_1D):
        _DataSet_double() 
        @staticmethod
        _DataSet * Alloc() 
        double& operator[](size_t idx)
        double& index_opr "operator[]"(size_t idx)
        #double operator[](size_t idx) const 
        const vector[double]& Data() const 
        #void operator =(const vector[double]& rhs)
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
        void WriteBuffer(_CpptrajFile&, size_t) const 
        void Append(const _DataSet_double&)
        #iterator begin() 
        #iterator end() 
        void SetNOE(double b, double bh, double r)
        double NOE_bound() const 
        double NOE_boundH() const 
        double NOE_rexp() const 
        void ShiftTorsions(double, double)


cdef class DataSet_double (DataSet_1D):
    cdef _DataSet_double* thisptr
    cdef bint py_free_mem 

