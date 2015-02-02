# distutils: language = c++
from DataSet_1D cimport *


cdef extern from "DataSet_Mesh.h": 
    cdef cppclass _DataSet_Mesh "DataSet_Mesh":
        _DataSet_Mesh() : _DataSet_1D(XYMESH, 12, 4)
        _DataSet_Mesh(int, double, double)
        _DataSet * Alloc() 
        size_t Size() const 
        int Sync() 
        void Info() const 
        int Allocate1D(size_t)
        void Add(size_t, const void *)
        double Dval(size_t idx) const 
        double Xcrd(size_t idx) const 
        void WriteBuffer(_CpptrajFile&, size_t) const 
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
