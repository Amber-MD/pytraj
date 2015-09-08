# distutils: language = c++
from libcpp.string cimport string
from ..Frame cimport _Frame, Frame
from ..core.cpptraj_core cimport _ArgList, ArgList
from ..Topology cimport _Topology, Topology
from ..cpptraj_dict cimport TrajFormatType


cdef extern from "TrajoutList.h": 
    cdef cppclass _Trajout "TrajoutList":
        _Trajout() 
        #int InitTrajWrite(const string&, _ArgList&, _Topology *, TrajFormatType)
        int InitTrajWrite "AddTrajout" (const string&, _ArgList&, _Topology *)
        void EndTraj "CloseTrajout"() 
        int WriteFrame "WriteTrajout"(int, const _Frame&)
        int SetupTrajWrite "SetupTrajout"(_Topology*)

cdef class Trajout:
    cdef _Trajout* thisptr
