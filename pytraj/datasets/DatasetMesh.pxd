# distutils: language = c++
from libcpp.vector cimport vector
from .dataset_1d cimport _DataSet_1D, DataSet_1D
from .base cimport _DataSet, DataSet


cdef extern from "DataSet_Mesh.h": 
    cdef cppclass _DataSet_Mesh "DataSet_Mesh" (_DataSet_1D):
        _DataSet_Mesh()
        _DataSet_Mesh(int, double, double)
        _DataSet * Alloc() 
        size_t Size() const 
        int Sync() 
        void Info() const 
        int Allocate1D(size_t)
        void Add(size_t, const void *)
        double Dval(size_t idx) const 
        double Xcrd(size_t idx) const 
        inline void AddXY(double, double)
        double X(int i) const 
        double Y(int i) const 
        void CalculateMeshX(int, double, double)
        int SetMeshXY(const _DataSet_1D&)
        double Integrate_Trapezoid(_DataSet_Mesh&) const 
        double Integrate_Trapezoid() const 
        int SetSplinedMeshY(const vector[double]&, const vector[double]&)
        int SetSplinedMesh(const _DataSet_1D&)
        int LinearRegression(double&, double&, double&, bint) const 

cdef class DataSet_Mesh(DataSet_1D):
    cdef _DataSet_Mesh* thisptr
    cdef public bint py_free_mem
