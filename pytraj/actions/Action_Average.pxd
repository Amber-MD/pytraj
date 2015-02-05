# distutils: language = c++
#from libcpp.string cimport string
from pytraj.actions.Action cimport _Action, Action, FunctPtr, _DispatchObject


cdef extern from "Action_Average.h": 
    cdef cppclass _Action_Average "Action_Average" (_Action):
        _Action_Average() 
        _DispatchObject * Alloc() 
        void Help() 
        #~_Action_Average() 


cdef class Action_Average (Action):
    cdef _Action_Average* thisptr

