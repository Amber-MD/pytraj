# distutils: language = c++
from pytraj.actions.Action cimport _Action, Action, FunctPtr, _DispatchObject
#from DataSet_Coords cimport *


cdef extern from "Action_CreateCrd.h": 
    cdef cppclass _Action_CreateCrd "Action_CreateCrd" (_Action):
        _Action_CreateCrd() 
        _DispatchObject * Alloc() 
        void Help() 


cdef class Action_CreateCrd (Action):
    cdef _Action_CreateCrd* thisptr

