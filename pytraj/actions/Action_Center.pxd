# distutils: language = c++
from pytraj.actions.Action cimport _Action, Action, FunctPtr, _DispatchObject


cdef extern from "Action_Center.h": 
    cdef cppclass _Action_Center "Action_Center" (_Action):
        _Action_Center() 
        _DispatchObject * Alloc() 
        void Help() 


cdef class Action_Center (Action):
    cdef _Action_Center* thisptr

