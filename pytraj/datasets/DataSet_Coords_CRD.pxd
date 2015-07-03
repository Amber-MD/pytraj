# distutils: language = c++

from .DataSet_Coords cimport _DataSet_Coords, DataSet_Coords
from .DataSet cimport _DataSet, DataSet
from .DataSet_1D cimport _DataSet_1D, DataSet_1D
from ..Frame cimport _Frame, Frame
from ..AtomMask cimport _AtomMask, AtomMask


cdef extern from "DataSet_Coords_CRD.h": 
    cdef cppclass _DataSet_Coords_CRD "DataSet_Coords_CRD" (_DataSet_Coords):
        _DataSet_Coords_CRD() 
        @staticmethod
        _DataSet * Alloc() 
        size_t Size() const 
        int Sync() 
        void Info() const 
        int Allocate1D(size_t)
        void Add(size_t, const void *)
        double Dval(size_t)const 
        double Xcrd(size_t idx)const 
        inline void AddFrame(const _Frame& fIn)
        inline void GetFrame(int idx, _Frame & fIn)
        inline void GetFrame(int idx, _Frame & fIn, const _AtomMask& mIn)
        inline void SetCRD(int idx, const _Frame& fIn)


cdef class DataSet_Coords_CRD (DataSet_Coords):
    cdef _DataSet_Coords_CRD* thisptr
