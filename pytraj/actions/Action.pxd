# distutils: language = c++
from libcpp.string cimport string
from pytraj.core.cpptraj_core cimport _DispatchObject, DispatchObject, FunctPtr
from pytraj.core.DataFileList cimport _DataFileList, DataFileList
from pytraj.ArgList cimport _ArgList, ArgList
from pytraj.datasets.DataSetList cimport _DataSetList, DataSetList
from pytraj.core.TopologyList cimport _TopologyList, TopologyList
from pytraj.Topology cimport _Topology, Topology
from pytraj.Frame cimport _Frame, Frame

cdef extern from "Action.h": 
    # Action.h
    ctypedef enum RetType "Action::RetType":
        OK "Action::OK"
        ERR "Action::ERR"
        USEORIGINALFRAME "Action::USEORIGINALFRAME"
        SUPPRESSCOORDOUTPUT "Action::SUPPRESSCOORDOUTPUT"
    cdef cppclass _Action "Action":
        #virtual ~_Action() 
        #RetType Init(_ArgList&, _TopologyList *, _FrameList *, _DataSetList *, _DataFileList *, int)
        RetType Init(_ArgList&, _TopologyList *, _DataSetList *, _DataFileList *, int)
        RetType Setup(_Topology *, _Topology * *)
        RetType DoAction(int, _Frame *, _Frame * *)
        void Print() 


cdef class Action:
    cdef _Action* baseptr
    cdef public int n_frames
