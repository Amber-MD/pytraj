# distutils: language = c++
from libcpp.string cimport string
from ..Topology cimport _Topology, Topology
from ..ArgList cimport _ArgList, ArgList
from ..Frame cimport _Frame, Frame
from ..AtomMask cimport _AtomMask, AtomMask
from ..FileName cimport _FileName, FileName
from .DataSet_Coords cimport _DataSet_Coords, DataSet_Coords
from .DataSet cimport _DataSet, DataSet
from .DataSet_1D cimport _DataSet_1D, DataSet_1D


cdef extern from "DataSet_Coords_REF.h": 
    cdef cppclass _DataSet_Coords_REF "DataSet_Coords_REF" (_DataSet_Coords):
        _DataSet_Coords_REF() 

        # turn off those methods since they are in parent class
        @staticmethod
        _DataSet * Alloc() 
        size_t Size() const 
        #int Sync() 
        #void Info() const 
        #void Add(size_t, const void *)
        #int AllocateCoords(size_t)
        #inline void AddFrame(const _Frame& fIn)
        #inline void GetFrame(int idx, _Frame& fIn)
        #inline void GetFrame(int idx, _Frame& fIn, const _AtomMask& mIn)
        #inline void SetCRD(int idx, const _Frame& fIn)

        int LoadRef(const string&, const _Topology&, int)
        int SetupRef_Frame(const string&, const string&, const _Topology&, _ArgList&, int)
        int SetupRef_Frame(_DataSet_Coords *, const string&, int, int)
        int StripRef(const string&)
        int StripRef(const _AtomMask&)
        const _Frame& RefFrame() const 
        const _FileName& _FrameFilename() const 
        int RefIndex() const 

cdef class DataSet_Coords_REF (DataSet_Coords):
    cdef _DataSet_Coords_REF* thisptr
