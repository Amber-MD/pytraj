# distutils: language = c++
from pytraj.actions.Action cimport _Action, Action, FunctPtr, _DispatchObject


cdef extern from "Action_FilterByData.h": 
    cdef cppclass _Action_FilterByData "Action_FilterByData" (_Action):
        _Action_FilterByData() 
        _DispatchObject * Alloc() 
        void Help() 


cdef class Action_FilterByData (Action):
    cdef _Action_FilterByData* thisptr

