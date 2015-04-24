# distutils: language = c++
from libcpp.vector cimport vector
from .DataSet cimport _DataSet, DataSet
from .DataSet_1D cimport _DataSet_1D, DataSet_1D
from ..math.Grid cimport _Grid
from ..math.Vec3 cimport _Vec3, Vec3
from ..CpptrajFile cimport _CpptrajFile, CpptrajFile


cdef extern from "DataSet_Vector.h": 
    cdef cppclass _DataSet_Vector "DataSet_Vector" (_DataSet_1D):
        _DataSet_Vector() 
        _DataSet * Alloc() 
        void WriteBuffer(_CpptrajFile&, size_t) const 
        void SetIred() 
        bint IsIred() const 
        void reset() 
        void Resize(size_t s)
        void Resize(size_t s, const _Vec3& v)
        bint Empty() const 
        #const _Vec3& operator[](int i) const 
        _Vec3& index_opr "operator[]" (int i)
        const _Vec3& OXYZ(int i) const 
        void ReserveVecs(size_t n)
        void AddVxyz(const _Vec3& v)
        void AddVxyz(const _Vec3& v, const _Vec3& c)
        #const_iterator begin() const 
        #const_iterator end() const 
        const _Vec3& Back() const 
        int CalcSphericalHarmonics(int)
        #const _ComplexArray& SphericalHarmonics(int) const 
        double SphericalHarmonicsNorm(int)


cdef class DataSet_Vector (DataSet_1D):
    cdef _DataSet_Vector* thisptr
    cdef bint py_free_mem

