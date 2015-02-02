# distutils: language = c++
#from libcpp.vector cimport vector
#from libcpp.string cimport string
from pytraj.actions.Action cimport _Action, Action, FunctPtr, _DispatchObject


cdef extern from "Action_Jcoupling.h": 
    cdef cppclass _Action_Jcoupling "Action_Jcoupling" (_Action):
        _Action_Jcoupling() 
        #~_Action_Jcoupling() 
        _DispatchObject * Alloc() 
        void Help() 


cdef class Action_Jcoupling (Action):
    cdef _Action_Jcoupling* thisptr
