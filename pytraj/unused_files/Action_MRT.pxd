# distutils: language = c++
#from libcpp.string cimport string
from pytraj.actions.Action cimport _Action, Action, FunctPtr, _DispatchObject


cdef extern from "Action_MRT.h": 
    cdef cppclass _Action_MRT "Action_MRT" (_Action):
        _Action_MRT() 
        _DispatchObject * Alloc() 
        void Help() 


cdef class Action_MRT (Action):
    cdef _Action_MRT* thisptr

