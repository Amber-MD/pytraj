# distutils: language = c++
from libcpp.vector cimport vector
from pytraj.datasets.DataSet cimport _DataSet, DataSet
from pytraj.datasets.DataSet_3D cimport _DataSet_3D, DataSet_3D
from pytraj.Grid cimport _Grid
from pytraj.Vec3 cimport _Vec3, Vec3
from pytraj.CpptrajFile cimport _CpptrajFile, CpptrajFile
from pytraj.ComplexArray cimport _ComplexArray, ComplexArray


cdef extern from "DataSet_Vector.h": 
    cdef cppclass _DataSet_Vector "DataSet_Vector":
        _DataSet_Vector() 
        _DataSet * Alloc() 
        size_t Size() const 
        int Sync() 
        void Info() const 
        int Allocate1D(size_t)
        inline void Add(size_t, const void *)
        double Dval(size_t) const 
        double Xcrd(size_t idx) const 
        void WriteBuffer(_CpptrajFile&, size_t) const 
        void SetIred() 
        bint IsIred() const 
        void reset() 
        void Resize(size_t s)
        void Resize(size_t s, const _Vec3& v)
        bint Empty() const 
        const _Vec3& operator[](int i) const 
        _Vec3& operator[](int i)
        const _Vec3& OXYZ(int i) const 
        void ReserveVecs(size_t n)
        void AddVxyz(const _Vec3& v)
        void AddVxyz(const _Vec3& v, const _Vec3& c)
        #const_iterator begin() const 
        #const_iterator end() const 
        const _Vec3& Back() const 
        int CalcSphericalHarmonics(int)
        const _ComplexArray& SphericalHarmonics(int) const 
        double SphericalHarmonicsNorm(int)


cdef class DataSet_Vector:
    cdef _DataSet_Vector* thisptr

