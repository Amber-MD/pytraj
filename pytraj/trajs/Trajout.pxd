# distutils: language = c++
from libcpp.string cimport string
from pytraj.trajs.TrajectoryFile cimport *
from pytraj.Range cimport *
#from pytraj.actions.ActionFrameCounter cimport *
from pytraj.Frame cimport _Frame, Frame


cdef extern from "Trajout.h": 
    cdef cppclass _Trajout "Trajout":
        _Trajout() 
        #~_Trajout() 
        inline int InitTrajWrite(const string&, _ArgList&, _Topology *, TrajFormatType)
        inline int InitTrajWrite(const string&, _Topology *, TrajFormatType)
        int InitStdoutTrajWrite(_ArgList&, _Topology *, _TrajFormatType)
        int InitEnsembleTrajWrite(const string&, const _ArgList&, _Topology *, TrajFormatType, int)
        int InitTrajWriteWithArgs(const string&, const char *, _Topology *, TrajFormatType)
        void EndTraj() 
        int WriteFrame(int, _Topology *, const _Frame&)
        void PrintInfo(int) const 
        bint TrajIsOpen() const 
        int NumFramesProcessed() const 


cdef class Trajout:
    cdef _Trajout* thisptr

