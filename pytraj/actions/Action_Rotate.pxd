# distutils: language = c++
from pytraj.actions.Action cimport _Action, Action, FunctPtr, _DispatchObject


cdef extern from "Action_Rotate.h": 
    cdef cppclass _Action_Rotate "Action_Rotate" (_Action):
        _Action_Rotate() 
        _DispatchObject * Alloc() 
        void Help() 


cdef class Action_Rotate (Action):
    cdef _Action_Rotate* thisptr

