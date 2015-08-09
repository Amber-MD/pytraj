# distutils: language = c++
from .DataSet cimport _DataSet, DataSet


cdef extern from "DataSet_1D.h": 
    cdef cppclass _DataSet_1D "DataSet_1D" (_DataSet):
        _DataSet_1D() 
        _DataSet_1D(_DataSet)
        # virtual methods
        #virtual ~_DataSet_1D() 
        int Allocate1D(size_t)
        double Dval(size_t) const
        double Xcrd(size_t) const
        # end virtual methods
        inline bint IsTorsionArray() const 
        double Avg() const 
        double Avg(double& sd) const 
        double Min() const 
        double Max() const 
        int CrossCorr(const _DataSet_1D&, _DataSet_1D&, int, bint, bint) const 
        double CorrCoeff(const _DataSet_1D&) const 


cdef class DataSet_1D (DataSet):
    # baseptr0 is from DataSet
    cdef _DataSet_1D* baseptr_1
