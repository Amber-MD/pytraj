# distutils: language = c++
from libcpp.vector cimport vector
from libcpp.string cimport string
from .DataSet cimport _DataSet, DataSet
from .DataSet_1D cimport _DataSet_1D, DataSet_1D
from .DataSet_Coords cimport _DataSet_Coords, DataSet_Coords
from ..trajs.Trajin cimport _Trajin, Trajin
from ..Frame cimport _Frame, Frame
from ..Topology cimport _Topology, Topology
from ..ArgList cimport _ArgList, ArgList
from ..AtomMask cimport _AtomMask, AtomMask


cdef extern from "DataSet_Coords_TRJ.h": 
    cdef cppclass _DataSet_Coords_TRJ "DataSet_Coords_TRJ" (_DataSet_Coords):
        _DataSet_Coords_TRJ() 
        #~_DataSet_Coords_TRJ() 
        @staticmethod
        _DataSet * Alloc() 
        int AddSingleTrajin(const string&, _ArgList&, _Topology *)
        int AddInputTraj(_Trajin *)
        size_t Size() const 
        int Allocate1D(size_t)
        void Add(size_t, const void *)
        double Dval(size_t) const 
        double Xcrd(size_t idx) const 
        void _AddFrame "AddFrame"(const _Frame& fIn)
        void SetCRD(int idx, const _Frame& fIn)
        void _GetFrame "GetFrame"(int idx, _Frame& fIn)
        void _GetFrame "GetFrame"(int idx, _Frame& fIn, const _AtomMask& mIn)
        int UpdateTrjFrames(int max_frames)


cdef class DataSet_Coords_TRJ(DataSet_Coords):
    cdef _DataSet_Coords_TRJ* thisptr
