# distutils: language = c++
from pytraj.actions.Action cimport _Action, Action, FunctPtr, _DispatchObject


cdef extern from "Action_NMRrst.h": 
    cdef cppclass _Action_NMRrst "Action_NMRrst" (_Action):
        _Action_NMRrst() 
        _DispatchObject * Alloc() 
        void Help() 


cdef class Action_NMRrst (Action):
    cdef _Action_NMRrst* thisptr
