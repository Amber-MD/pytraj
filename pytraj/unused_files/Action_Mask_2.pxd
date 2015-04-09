# distutils: language = c++
#from libcpp.string cimport string
from pytraj.actions.Action cimport _Action, Action, FunctPtr, _DispatchObject
#from TrajectoryFile cimport *


cdef extern from "src/Action_Mask_2_cpp.h": 
    cdef cppclass _Action_Mask_2 "Action_Mask_2" (_Action):
        _Action_Mask_2() 
        _DispatchObject * Alloc() 
        void Help() 


cdef class Action_Mask_2 (Action):
    cdef _Action_Mask_2* thisptr

