# distutils: language = c++
#from libcpp.vector cimport vector
#from libcpp.string cimport string
from pytraj.actions.Action cimport _Action, Action, FunctPtr, _DispatchObject
#from AxisType cimport *
#from Range cimport *
#from DataSet_1D cimport *


cdef extern from "Action_NAstruct.h": 
    cdef cppclass _Action_NAstruct "Action_NAstruct" (_Action):
        _Action_NAstruct() 
        #~_Action_NAstruct() 
        _DispatchObject * Alloc() 
        void Help() 


cdef class Action_NAstruct (Action):
    cdef _Action_NAstruct* thisptr

