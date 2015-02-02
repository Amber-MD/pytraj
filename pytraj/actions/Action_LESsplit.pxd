# distutils: language = c++
#from libcpp.vector cimport vector
#from libcpp.string cimport string
from pytraj.actions.Action cimport _Action, Action, FunctPtr, _DispatchObject
#from Trajout cimport *


cdef extern from "Action_LESsplit.h": 
    cdef cppclass _Action_LESsplit "Action_LESsplit" (_Action):
        _Action_LESsplit() 
        #~_Action_LESsplit() 
        _DispatchObject * Alloc() 
        void Help() 


cdef class Action_LESsplit (Action):
    cdef _Action_LESsplit* thisptr

