# distutils: language = c++
#from libcpp.vector cimport vector
#from libcpp.string cimport string
from pytraj.actions.Action cimport _Action, Action, FunctPtr, _DispatchObject
#from ImagedAction cimport *


cdef extern from "Action_Watershell.h": 
    cdef cppclass _Action_Watershell "Action_Watershell" (_Action):
        _Action_Watershell() 
        #~_Action_Watershell() 
        _DispatchObject * Alloc() 
        void Help() 


cdef class Action_Watershell (Action):
    cdef _Action_Watershell* thisptr

