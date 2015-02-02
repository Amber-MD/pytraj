# distutils: language = c++
#from libcpp.vector cimport vector
from pytraj.actions.Action cimport _Action, Action, FunctPtr, _DispatchObject
#from ImagedAction cimport *


cdef extern from "Action_LIE.h": 
    cdef cppclass _Action_LIE "Action_LIE" (_Action):
        _Action_LIE() 
        _DispatchObject * Alloc() 
        void Help() 


cdef class Action_LIE (Action):
    cdef _Action_LIE* thisptr

