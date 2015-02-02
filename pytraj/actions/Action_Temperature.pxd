# distutils: language = c++
from pytraj.actions.Action cimport _Action, Action, FunctPtr, _DispatchObject


cdef extern from "Action_Temperature.h": 
    cdef cppclass _Action_Temperature "Action_Temperature" (_Action):
        _Action_Temperature() 
        _DispatchObject * Alloc() 
        void Help() 


cdef class Action_Temperature (Action):
    cdef _Action_Temperature* thisptr

