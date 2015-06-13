# distutils: language = c++
from __future__ import absolute_import
from .DataSet cimport _DataSet, DataSet
from .DataSet_3D cimport _DataSet_3D, DataSet_3D
from ..math.Grid cimport _Grid
from ..math.Vec3 cimport _Vec3, Vec3


cdef extern from "DataSet_GridFlt.h": 
    cdef cppclass _DatasetGridFloat "DataSet_GridFlt" (_DataSet_3D):
        _DatasetGridFloat()
        float& index_opr "operator[]"(size_t idx)
        _DataSet * Alloc() 
        const _Grid[float]& InternalGrid() const 
        size_t Size() const 
        int Sync() 
        void Info() const 
        int Allocate3D(size_t x, size_t y, size_t z)
        double GetElement(int x, int y, int z) const 
        void SetElement(int x, int y, int z, float v)
        double operator[](size_t idx) const 
        size_t NX() const 
        size_t NY() const 
        size_t NZ() const 
        #iterator begin() 
        #iterator end() 
        inline long int Increment(const _Vec3&, float)
        inline long int Increment(const double *, float)
        inline long int Increment(int, int, int, float)
        float GridVal(int x, int y, int z) const 
        long int CalcIndex(int i, int j, int k) const 

cdef class DatasetGridFloat (DataSet_3D):
    cdef _DatasetGridFloat* thisptr
    cdef public bint py_free_mem
