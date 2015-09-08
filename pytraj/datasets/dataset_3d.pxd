# distutils: language = c++
from .base cimport _DataSet, DataSet, DataType
from ..math.Vec3 cimport _Vec3, Vec3
from .base cimport _DataSet, DataSet
from ..math.Grid cimport _Grid
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
