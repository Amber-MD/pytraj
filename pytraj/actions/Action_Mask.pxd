# distutils: language = c++
#from libcpp.string cimport string
from pytraj.actions.Action cimport _Action, Action, FunctPtr, _DispatchObject
#from TrajectoryFile cimport *


cdef extern from "Action_Mask.h": 
    cdef cppclass _Action_Mask "Action_Mask" (_Action):
        _Action_Mask() 
        _DispatchObject * Alloc() 
        void Help() 


cdef class Action_Mask (Action):
    cdef _Action_Mask* thisptr

