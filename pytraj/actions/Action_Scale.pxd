# distutils: language = c++
from pytraj.actions.Action cimport _Action, Action, FunctPtr, _DispatchObject


cdef extern from "Action_Scale.h": 
    cdef cppclass _Action_Scale "Action_Scale" (_Action):
        _Action_Scale() 
        _DispatchObject * Alloc() 
        void Help() 


cdef class Action_Scale (Action):
    cdef _Action_Scale* thisptr

