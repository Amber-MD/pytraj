# distutils: language = c++
from libcpp.string cimport string
from ..frame cimport _Frame, Frame
from ..core.c_core cimport _ArgList, ArgList, Box
from ..topology cimport _Topology, Topology
from ..c_dict cimport TrajFormatType
from ..core.coordinfo cimport _CoordinateInfo, CoordinateInfo


cdef extern from "TrajoutList.h": 
    cdef cppclass _Trajout "TrajoutList":
        _Trajout() 
        #int InitTrajWrite(const string&, _ArgList&, _Topology *, TrajFormatType)
        int InitTrajWrite "AddTrajout" (const string&, _ArgList&, _Topology *)
        void EndTraj "CloseTrajout"() 
        int WriteFrame "WriteTrajout"(int, const _Frame&)
        int SetupTrajWrite "SetupTrajout"(_Topology*, _CoordinateInfo, int)

cdef class TrajectoryWriter:
    cdef _Trajout* thisptr
    cdef unsigned int count
