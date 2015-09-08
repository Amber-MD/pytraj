# distutils: language = c++
from libcpp.vector cimport vector
from libcpp.string cimport string
from ..Frame cimport _Frame, Frame
from ..Topology cimport _Topology, Topology
from ..core.cpptraj_core cimport _AtomMask, AtomMask, _ArgList, ArgList


cdef extern from "DataSet_Coords_TRJ.h" nogil: 
    cdef cppclass _TrajectoryCpptraj "DataSet_Coords_TRJ":
        _TrajectoryCpptraj() 
        int AddSingleTrajin(const string&, _ArgList&, _Topology *)
        size_t Size() const 
        void GetFrame(int idx, _Frame& fIn)
        void GetFrame(int idx, _Frame& fIn, _AtomMask& mIn)
        void SetTopology(const _Topology&)
        inline const _Topology& Top() const 


cdef class TrajectoryCpptraj:
    cdef Topology _top
    cdef _TrajectoryCpptraj* thisptr
    cdef object tmpfarray
    cdef list _filelist
