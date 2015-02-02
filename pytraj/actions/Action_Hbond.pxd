# distutils: language = c++
#from libcpp.vector cimport vector
#from libcpp.string cimport string
from pytraj.actions.Action cimport _Action, Action, FunctPtr, _DispatchObject
#from ImagedAction cimport *
#from DataSet_integer cimport *


cdef extern from "Action_Hbond.h": 
    cdef cppclass _Action_Hbond "Action_Hbond" (_Action):
        _Action_Hbond() 
        _DispatchObject * Alloc() 
        void Help() 


cdef class Action_Hbond (Action):
    cdef _Action_Hbond* thisptr

