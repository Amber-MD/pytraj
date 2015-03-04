# distutils: language = c++
from libcpp.string cimport string
from pytraj.trajs.TrajectoryFile cimport *
from pytraj.Range cimport *
from pytraj.Frame cimport _Frame, Frame


cdef extern from "Trajout_Single.h": 
    cdef cppclass _Trajout "Trajout_Single":
        _Trajout() 
        #~_Trajout() 
        int InitTrajWrite(const string&, _ArgList&, _Topology *, TrajFormatType)
        void EndTraj() 
        int WriteFrame "WriteSingle"(int, const _Frame&)
        #int WriteFrame(int, _Topology *, const _Frame&)

cdef class Trajout:
    cdef _Trajout* thisptr
