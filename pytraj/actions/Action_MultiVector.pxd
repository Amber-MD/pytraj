# distutils: language = c++
from pytraj.actions.Action cimport _Action, Action, FunctPtr, _DispatchObject


cdef extern from "Action_MultiVector.h": 
    cdef cppclass _Action_MultiVector "Action_MultiVector" (_Action):
        _Action_MultiVector() 
        _DispatchObject * Alloc() 
        void Help() 


cdef class Action_MultiVector (Action):
    cdef _Action_MultiVector* thisptr

