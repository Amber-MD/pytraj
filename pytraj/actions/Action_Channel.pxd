# distutils: language = c++
from pytraj.actions.Action cimport _Action, Action, FunctPtr, _DispatchObject


cdef extern from "Action_Channel.h": 
    cdef cppclass _Action_Channel "Action_Channel" (_Action):
        _Action_Channel() 
        _DispatchObject * Alloc() 
        void Help() 


cdef class Action_Channel (Action):
    cdef _Action_Channel* thisptr

