# distutils: language = c++
#from libcpp.vector cimport vector
#from libcpp.string cimport string
from pytraj.actions.Action cimport _Action, Action, FunctPtr, _DispatchObject
#from ImagedAction cimport *
#from OnlineVarT cimport *


cdef extern from "Action_OrderParameter.h": 
    cdef cppclass _Action_OrderParameter "Action_OrderParameter" (_Action):
        _Action_OrderParameter() 
        _DispatchObject * Alloc() 
        void Help() 


cdef class Action_OrderParameter (Action):
    cdef _Action_OrderParameter* thisptr

