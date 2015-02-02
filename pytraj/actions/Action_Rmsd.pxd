# distutils: language = c++
#from pytraj.actions.Action cimport _Action, Action
from pytraj.actions.Action cimport _Action, Action, FunctPtr, _DispatchObject

cdef extern from "Action_Rmsd.h": 
    cdef cppclass _Action_Rmsd"Action_Rmsd" (_Action):
        _Action_Rmsd() 
        _DispatchObject * Alloc() 
        void Help() 

cdef class Action_Rmsd(Action):
    cdef _Action_Rmsd* thisptr
    #cdef DispatchAllocatorType func
