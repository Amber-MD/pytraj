# distutils: language = c++
#from libcpp.vector cimport vector
#from pytraj.actions.Action cimport *
from pytraj.actions.Action cimport _Action, Action, FunctPtr, _DispatchObject
#from Array1D cimport *


cdef extern from "Action_FilterByData.h": 
    cdef cppclass _Action_FilterByData "Action_FilterByData" (_Action):
        _Action_FilterByData() 
        _DispatchObject * Alloc() 
        void Help() 
        size_t Determine_Frames() const 
        #RetType Init(_ArgList&, _TopologyList *, _FrameList *, _DataSetList *, _DataFileList *, int)
        #RetType DoAction(int, _Frame *, _Frame * *)


cdef class Action_FilterByData (Action):
    cdef _Action_FilterByData* thisptr

