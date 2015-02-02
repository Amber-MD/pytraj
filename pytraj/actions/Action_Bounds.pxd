# distutils: language = c++
#from libcpp.string cimport string
from pytraj.actions.Action cimport _Action, Action, FunctPtr, _DispatchObject


cdef extern from "Action_Bounds.h": 
    cdef cppclass _Action_Bounds "Action_Bounds" (_Action):
        _Action_Bounds() 
        _DispatchObject * Alloc() 
        void Help() 


cdef class Action_Bounds (Action):
    cdef _Action_Bounds* thisptr

