# distutils: language = c++
from pytraj.actions.Action cimport _Action, Action, FunctPtr, _DispatchObject


cdef extern from "Action_Energy.h": 
    cdef cppclass _Action_Energy "Action_Energy" (_Action):
        _Action_Energy() 
        _DispatchObject * Alloc() 
        void Help() 


cdef class Action_Energy (Action):
    cdef _Action_Energy* thisptr

