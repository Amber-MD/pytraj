# distutils: language = c++
#from libcpp.vector cimport vector
from pytraj.actions.Action cimport _Action, Action, FunctPtr, _DispatchObject


cdef extern from "Action_Pucker.h": 
    cdef cppclass _Action_Pucker "Action_Pucker" (_Action):
        _Action_Pucker() 
        _DispatchObject * Alloc() 
        void Help() 


cdef class Action_Pucker (Action):
    cdef _Action_Pucker* thisptr

