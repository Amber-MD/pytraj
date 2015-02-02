# distutils: language = c++
#from libcpp.vector cimport vector
from pytraj.actions.Action cimport _Action, Action, FunctPtr, _DispatchObject


cdef extern from "Action_Surf.h": 
    cdef cppclass _Action_Surf "Action_Surf" (_Action):
        _Action_Surf() 
        _DispatchObject * Alloc() 
        void Help() 


cdef class Action_Surf (Action):
    cdef _Action_Surf* thisptr

