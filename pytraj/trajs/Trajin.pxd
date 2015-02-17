# distutils: language = c++
#from pytraj.trajs.TrajectoryFile cimport *
from libcpp.string cimport string
from pytraj.ArgList cimport _ArgList, ArgList
from pytraj.Frame cimport _Frame, Frame
from pytraj.Box cimport _Box, Box
from pytraj.trajs.TrajectoryIO cimport _TrajectoryIO, TrajectoryIO
from pytraj.trajs.TrajectoryFile cimport _TrajectoryFile, TrajectoryFile
from pytraj.Topology cimport _Topology, Topology


cdef extern from "Trajin.h": 
    cdef cppclass _Trajin "Trajin" (_TrajectoryFile):
        _Trajin()        

        #virtual ~_Trajin() 
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
        int SetupTrajIO(const string&, _TrajectoryIO &, _ArgList &)
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

cdef class Trajin (TrajectoryFile):
    #( baseptr0 is from TrajectoryFile
    # create tmpfarray to hold sub FrameArray
    # traj[0:10][0] will give wrong answer
    cdef object tmpfarray
    cdef _Trajin* baseptr_1
    cdef public bint debug
