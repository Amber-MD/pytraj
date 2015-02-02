# distutils: language = c++
#from libcpp.vector cimport vector
from pytraj.actions.Action cimport _Action, Action, FunctPtr, _DispatchObject
from pytraj.actions.Action cimport _ArgList, _TopologyList, RetType
from pytraj.actions.Action cimport _FrameList, _DataSetList, _DataFileList
from pytraj.Frame cimport _Frame, Frame
from pytraj.Topology cimport _Topology, Topology
#from ImagedAction cimport *


cdef extern from "Action_CheckStructure.h": 
#    cdef cppclass _Action_CheckStructure::Problem "Action_CheckStructure::Problem" (_Action):
#        Problem() 
#        Problem(int f)
#        Problem(const Problem& rhs)


    cdef cppclass _Action_CheckStructure "Action_CheckStructure" (_Action):
        _Action_CheckStructure() 
        _DispatchObject * Alloc() 
        void Help() 
        #~_Action_CheckStructure() 
        RetType Init(_ArgList&, _TopologyList *, _FrameList *, _DataSetList *, _DataFileList *, int)
        RetType Setup(_Topology *, _Topology **)
        int Check_Frame(int, const _Frame&)
        void Print() 


cdef class Action_CheckStructure (Action):
    cdef _Action_CheckStructure* thisptr
