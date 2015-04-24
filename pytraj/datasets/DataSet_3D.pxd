# distutils: language = c++
from .DataSet cimport _DataSet, DataSet, DataType
from ..CpptrajFile cimport _CpptrajFile, CpptrajFile
from ..math.Vec3 cimport _Vec3, Vec3


cdef extern from "DataSet_3D.h": 
    cdef cppclass _DataSet_3D "DataSet_3D" (_DataSet):
        _DataSet_3D() 
        _DataSet_3D(DataType tIn, int wIn, int pIn)
        void Add(size_t, const void *)
        int Allocate_N_O_D(size_t, size_t, size_t, const _Vec3&, const _Vec3&)
        int Allocate_N_C_D(size_t, size_t, size_t, const _Vec3&, const _Vec3&)
        int Allocate_X_C_D(const _Vec3&, const _Vec3&, const _Vec3&)
        inline bint CalcBins(double, double, double, int&, int&, int&) const 
        inline double DX() const 
        inline double DY() const 
        inline double DZ() const 
        inline double OX() const 
        inline double OY() const 
        inline double OZ() const 
        inline double MX() const 
        inline double MY() const 
        inline double MZ() const 
        inline _Vec3 BinCorner(int, int, int)
        inline _Vec3 BinCenter(int, int, int)

cdef class DataSet_3D (DataSet):
    cdef _DataSet_3D* baseptr_1
