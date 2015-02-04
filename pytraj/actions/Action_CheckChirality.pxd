# distutils: language = c++
#from pytraj.actions.Action cimport *
from pytraj.actions.Action cimport _Action, Action, FunctPtr, _DispatchObject

cdef extern from "Action_CheckChirality.h": 
    cdef cppclass _Action_CheckChirality "Action_CheckChirality" (_Action):
        _Action_CheckChirality() 
        _DispatchObject * Alloc() 
        void Help() 


cdef class Action_CheckChirality (Action):
    cdef _Action_CheckChirality* thisptr

