# distutils: language = c++
from pytraj.actions.Action cimport _Action, Action, FunctPtr, _DispatchObject


cdef extern from "Action_Translate.h": 
    cdef cppclass _Action_Translate "Action_Translate" (_Action):
        _Action_Translate() 
        _DispatchObject * Alloc() 
        void Help() 


cdef class Action_Translate (Action):
    cdef _Action_Translate* thisptr

