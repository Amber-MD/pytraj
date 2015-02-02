# distutils: language = c++
from libcpp.string cimport string
from pytraj.Topology cimport _Topology, Topology
from pytraj.ReplicaDimArray cimport * 
from pytraj.CpptrajFile cimport _CpptrajFile, CpptrajFile
from pytraj.ArgList cimport _ArgList, ArgList
from pytraj.BaseIOtype cimport _BaseIOtype, BaseIOtype
from pytraj.Box cimport _Box, Box
from pytraj.Frame cimport _Frame, Frame


cdef extern from "TrajectoryIO.h": 
    cdef cppclass _TrajectoryIO "TrajectoryIO" (_BaseIOtype):
        _TrajectoryIO() 
        #virtual ~_TrajectoryIO() 

        # virtual methods
        bint ID_TrajFormat(_CpptrajFile&) 
        int setupTrajin(const string&, _Topology *) 
        int setupTrajout(const string&, _Topology *, int, bint) 
        int openTrajin()  
        int readFrame(int, _Frame&) 
        int readVelocity(int, _Frame&) 
        int writeFrame(int, const _Frame&) 
        void closeTraj()  
        void Info()  
        int processWriteArgs(_ArgList&) 
        int processReadArgs(_ArgList&) 
        # end virtual methods

        bint HasBox() const 
        const _Box& TrajBox() const 
        bint HasV() const 
        bint HasT() const 
        const string& Title() const 
        const _ReplicaDimArray& ReplicaDimensions() const 
        void SetDebug(int dIn)
        void SetBox(const _Box& bIn)
        void SetVelocity(bint vIn)
        void SetTemperature(bint tIn)
        void SetTitle(const string& tIn)
        void SetReplicaDims(const _ReplicaDimArray& rIn)


cdef class TrajectoryIO (BaseIOtype):
    cdef _TrajectoryIO* baseptr_1
