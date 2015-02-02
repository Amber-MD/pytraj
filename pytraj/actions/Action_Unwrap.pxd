# distutils: language = c++
from pytraj.actions.Action cimport _Action, Action, FunctPtr, _DispatchObject
#from ImageTypes cimport *


cdef extern from "Action_Unwrap.h": 
    cdef cppclass _Action_Unwrap "Action_Unwrap" (_Action):
        _Action_Unwrap() 
        _DispatchObject * Alloc() 
        void Help() 


cdef class Action_Unwrap (Action):
    cdef _Action_Unwrap* thisptr

