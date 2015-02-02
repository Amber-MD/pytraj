# distutils: language = c++
#from libcpp.vector cimport vector
#from libcpp.string cimport string
from pytraj.actions.Action cimport _Action, Action, FunctPtr, _DispatchObject


cdef extern from "Action_FixAtomOrder.h": 
    cdef cppclass _Action_FixAtomOrder "Action_FixAtomOrder" (_Action):
        _Action_Fix_AtomOrder() 
        #~_Action_Fix_AtomOrder() 
        _DispatchObject * Alloc() 
        void Help() 


cdef class Action_FixAtomOrder (Action):
    cdef _Action_FixAtomOrder* thisptr

