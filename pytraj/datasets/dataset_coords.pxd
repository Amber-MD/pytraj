# distutils: language = c++
from libcpp.string  cimport string
from .base cimport _DataSet, DataSet
from ..Frame cimport _Frame, Frame
from ..Topology cimport _Topology, Topology
from ..core.cpptraj_core cimport _ArgList, ArgList, _AtomMask, AtomMask


cdef extern from "DataSet_Coords.h": 
    cdef cppclass _DataSet_Coords "DataSet_Coords" (_DataSet):
        _DataSet_Coords() 
        _DataSet_Coords(DataType)
        #virtual ~_DataSet_Coords() 
        _Frame AllocateFrame() const 
        
        # virtual methods
        void AddFrame(const _Frame&) 
        void SetCRD(int, const _Frame&) 
        void GetFrame(int, _Frame&) 
        void GetFrame(int, _Frame&, const _AtomMask&) 
        # end virtual methods

        void SetTopology(const _Topology&)
        inline const _Topology& Top() const 


cdef class DataSet_Coords (DataSet):
    # DataSet has baseptr0
    cdef _DataSet_Coords* baseptr_1
    cdef Topology _top
    cdef bint py_free_mem

    # use tmpfarray object to hold Frame or Trajectory 
    # (if we want to use dset[0][0] correctly)
    cdef object tmpfarray
# distutils: language = c++

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
        int RefIndex() const 

cdef class DataSet_Coords_REF (DataSet_Coords):
    cdef _DataSet_Coords_REF* thisptr
