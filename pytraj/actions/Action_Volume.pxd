# distutils: language = c++
from pytraj.actions.Action cimport _Action, Action, FunctPtr, _DispatchObject


cdef extern from "Action_Volume.h": 
    cdef cppclass _Action_Volume "Action_Volume" (_Action):
        _Action_Volume() 
        _DispatchObject * Alloc() 
        void Help() 


cdef class Action_Volume (Action):
    cdef _Action_Volume* thisptr

