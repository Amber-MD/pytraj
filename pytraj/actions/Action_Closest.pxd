# distutils: language = c++
#from libcpp.vector cimport vector
#from libcpp.string cimport string
from pytraj.actions.Action cimport _Action, Action, FunctPtr, _DispatchObject
#from ImagedAction cimport *


cdef extern from "Action_Closest.h": 
    cdef cppclass _Action_Closest "Action_Closest" (_Action):
        _Action_Closest() 
        _DispatchObject * Alloc() 
        void Help() 
        #~_Action_Closest() 


cdef class Action_Closest (Action):
    cdef _Action_Closest* thisptr
