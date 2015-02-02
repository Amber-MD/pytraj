# distutils: language = c++
#from libcpp.string cimport string
#from pytraj.actions.Action cimport *
from pytraj.actions.Action cimport _Action, Action, FunctPtr, _DispatchObject


cdef extern from "Action_Strip.h": 
    cdef cppclass _Action_Strip "Action_Strip" (_Action):
        _Action_Strip() 
        _DispatchObject * Alloc() 
        void Help() 
        #~_Action_Strip() 


    cdef cppclass _Action_Unstrip "Action_Unstrip" (_Action):
        _Action_Unstrip() 
        _DispatchObject * Alloc() 
        void Help() 


cdef class Action_Strip (Action):
    cdef _Action_Strip* thisptr

cdef class Action_Unstrip (Action):
    cdef _Action_Unstrip* thisptr

