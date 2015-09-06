# distutils: language = c++
from libcpp.vector cimport vector
from .base cimport DataSet, _DataSet
from .dataset_2d cimport _DataSet_2D, DataSet_2D
from ..Frame cimport *
from ..analyses.Analysis cimport *


cdef extern from "DataSet_Modes.h": 
    cdef cppclass _DataSet_Modes "DataSet_Modes" (_DataSet):
        _DataSet_Modes() 
        _DataSet * Alloc() 
        size_t Size() const 
        int Sync() 
        void Info() const 
        void Add(size_t, const void *)
        const Darray& AvgCrd() const 
        const Darray& Mass() const 
        int NavgCrd() const 
        double * Avg_FramePtr() 
        const double * Avg_FramePtr() const 
        void AllocateAvgCoords(int n)
        void SetAvgCoords(const _DataSet_2D&)
        int SetModes(bint, int, int, const double *, const double *)
        int CalcEigen(const _DataSet_2D&, int)
        void PrintModes() 
        int EigvalToFreq(double)
        int MassWtEigvect(Darray&)
        int ReduceVectors() 
        double Eigenvalue(int i) const 
        const double * Eigenvectors() const 
        const double * Eigenvector(int i) const 
        int Nmodes() const 
        int VectorSize() const 
        #MatrixType Type() const 
        bint IsReduced() const 

cdef class DataSet_Modes (DataSet):
    cdef _DataSet_Modes* thisptr
