# distutils: language = c++
from pytraj.actions.Action cimport _Action, Action, FunctPtr, _DispatchObject


cdef extern from "Action_CheckStructure.h": 
    cdef cppclass _Action_CheckStructure "Action_CheckStructure" (_Action):
        _Action_CheckStructure() 
        _DispatchObject * Alloc() 
        void Help() 


cdef class Action_CheckStructure (Action):
    cdef _Action_CheckStructure* thisptr
