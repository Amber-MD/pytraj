# distutils: language = c++
from libcpp.string cimport string
from pytraj.DispatchObject cimport _DispatchObject, DispatchObject
from pytraj.ArgList cimport _ArgList, ArgList
from pytraj.DataFileList cimport _DataFileList, DataFileList
from pytraj.DataSetList cimport _DataSetList, DataSetList
from pytraj.FrameList cimport _FrameList, FrameList
from pytraj.TopologyList cimport _TopologyList, TopologyList
from pytraj._FunctPtr cimport FunctPtr
from pytraj.Topology cimport _Topology, Topology
from pytraj.Frame cimport _Frame, Frame
from pytraj.FrameArray cimport FrameArray

cdef extern from "Action.h": 
    # Action.h
    ctypedef enum RetType "Action::RetType":
        OK "Action::OK"
        ERR "Action::ERR"
        USEORIGINALFRAME "Action::USEORIGINALFRAME"
        SUPPRESSCOORDOUTPUT "Action::SUPPRESSCOORDOUTPUT"
    cdef cppclass _Action "Action":
        #virtual ~_Action() 
        RetType Init(_ArgList&, _TopologyList *, _FrameList *, _DataSetList *, _DataFileList *, int)
        RetType Setup(_Topology *, _Topology * *)
        RetType DoAction(int, _Frame *, _Frame * *)
        void Print() 


cdef class Action:
    cdef _Action* baseptr
