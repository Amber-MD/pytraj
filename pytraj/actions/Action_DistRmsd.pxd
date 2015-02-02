# distutils: language = c++
from pytraj.actions.Action cimport _Action, Action, FunctPtr, _DispatchObject
#from ReferenceAction cimport *


cdef extern from "Action_DistRmsd.h": 
    cdef cppclass _Action_DistRmsd "Action_DistRmsd" (_Action):
        _Action_DistRmsd() 
        _DispatchObject * Alloc() 
        void Help() 


cdef class Action_DistRmsd (Action):
    cdef _Action_DistRmsd* thisptr

