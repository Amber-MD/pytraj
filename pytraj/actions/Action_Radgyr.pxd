# distutils: language = c++
from pytraj.actions.Action cimport _Action, Action, FunctPtr, _DispatchObject


cdef extern from "Action_Radgyr.h": 
    cdef cppclass _Action_Radgyr "Action_Radgyr" (_Action):
        _Action_Radgyr() 
        _DispatchObject * Alloc() 
        void Help() 


cdef class Action_Radgyr (Action):
    cdef _Action_Radgyr* thisptr

