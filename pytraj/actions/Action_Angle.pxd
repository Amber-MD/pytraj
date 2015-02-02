# distutils: language = c++
from pytraj.actions.Action cimport _Action, Action, FunctPtr, _DispatchObject


cdef extern from "Action_Angle.h": 
    cdef cppclass _Action_Angle "Action_Angle" (_Action):
        _Action_Angle() 
        _DispatchObject * Alloc() 
        void Help() 


cdef class Action_Angle (Action):
    cdef _Action_Angle* thisptr

