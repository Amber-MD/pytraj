# distutils: language = c++
from libcpp.string cimport string
from pytraj.ArgList cimport _ArgList, ArgList
from pytraj.Frame cimport _Frame, Frame
from pytraj.Topology cimport _Topology, Topology
from pytraj.core.FileName cimport _FileName, FileName
from pytraj.core.Box cimport _Box, Box
from pytraj.core.CoordinateInfo cimport _CoordinateInfo, CoordinateInfo

cdef extern from "Trajin.h": 
    cdef cppclass _Trajin "Trajin" nogil:
        _Trajin()        

        # star virtual methods
        int SetupTrajRead(const string&, _ArgList&, _Topology *)
        int ReadTrajFrame(int, _Frame&)
        int BeginTraj(bint)
        void EndTraj() 
        void PrintInfo(int) const 
        bint HasVelocity() const 
        int NreplicaDimension() const 
        # end virtual methods

        @staticmethod
        int CheckFrameArgs(_ArgList &, int, int &, int &, int &)
        inline bint CheckFinished() 
        inline void UpdateCounters() 
        inline int GetNextFrame(_Frame &)
        inline void SetTotalFrames(int)
        int CheckBoxInfo(const char *, _Box &, const _Box&)const 
        int setupFrameInfo() 
        void PrepareForRead(bint)
        void PrintInfoLine() const 
        void PrintFrameInfo() const 
        int TotalFrames() const 
        int TotalReadFrames() const 
        int CurrentFrame() const 
        int Start() const 
        int Offset() const 
        int NumFramesProcessed() const 
        bint IsEnsemble() const 
        void SetEnsemble(bint b)
        _CoordinateInfo TrajCoordInfo()

        void SetTrajFileName(const string&, bint)
        int SetTrajParm(_Topology *)
        _Topology * TrajParm ()const 
        const _FileName & TrajFilename ()const 

cdef class Trajin:
    cdef object tmpfarray
    cdef _Trajin* baseptr_1
    cdef public bint debug
    cdef Topology _top
    cdef public object _tmpobj
