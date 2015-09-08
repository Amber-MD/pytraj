# distutils: language = c++
from libcpp.string cimport string
from ..core.cpptraj_core cimport (_DispatchObject, DispatchObject,  FunctPtr)
from ..datafiles.datafiles cimport  _DataFileList, DataFileList
from ..core.TopologyList cimport _TopologyList, TopologyList
from ..core.cpptraj_core cimport _ArgList, ArgList
from ..datasets.DataSetList cimport _DataSetList, DataSetList
from ..Topology cimport _Topology, Topology
from ..Frame cimport _Frame, Frame


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
