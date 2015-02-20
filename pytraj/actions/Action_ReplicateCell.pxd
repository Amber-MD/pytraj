# distutils: language = c++
from pytraj.actions.Action cimport _Action, Action, FunctPtr, _DispatchObject


cdef extern from "Action_ReplicateCell.h": 
    cdef cppclass _Action_ReplicateCell "Action_ReplicateCell" (_Action):
        _Action_ReplicateCell() 
        _DispatchObject * Alloc() 
        void Help() 


cdef class Action_ReplicateCell (Action):
    cdef _Action_ReplicateCell* thisptr

