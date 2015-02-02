# distutils: language = c++
#from libcpp.vector cimport vector
#from libcpp.string cimport string
from pytraj.actions.Action cimport _Action, Action, FunctPtr, _DispatchObject
#from ImagedAction cimport *
#from OnlineVarT cimport *


cdef extern from "Action_Density.h": 
    cdef cppclass _Action_Density "Action_Density" (_Action):
        _Action_Density() 
        _DispatchObject * Alloc() 
        void Help() 


cdef class Action_Density (Action):
    cdef _Action_Density* thisptr

