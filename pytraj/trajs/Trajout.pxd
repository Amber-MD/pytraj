# distutils: language = c++
from libcpp.string cimport string
from pytraj.Frame cimport _Frame, Frame
from pytraj.ArgList cimport _ArgList, ArgList
from pytraj.Topology cimport _Topology, Topology
from pytraj.cpptraj_dict cimport TrajFormatType


cdef extern from "Trajout_Single.h": 
    cdef cppclass _Trajout "Trajout_Single":
        _Trajout() 
        int InitTrajWrite(const string&, _ArgList&, _Topology *, TrajFormatType)
        void EndTraj() 
        int WriteFrame "WriteSingle"(int, const _Frame&)
        int SetupTrajWrite(_Topology*)

cdef class Trajout:
    cdef _Trajout* thisptr
