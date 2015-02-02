# distutils: language = c++
from pytraj.actions.Action cimport _Action, Action, FunctPtr, _DispatchObject
#from ImagedAction cimport *


cdef extern from "Action_Radial.h": 
    cdef cppclass _Action_Radial "Action_Radial" (_Action):
        _Action_Radial() 
        #~_Action_Radial() 
        _DispatchObject * Alloc() 
        void Help() 


cdef class Action_Radial (Action):
    cdef _Action_Radial* thisptr

