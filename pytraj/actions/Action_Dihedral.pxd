# distutils: language = c++
from pytraj.actions.Action cimport _Action, Action, FunctPtr, _DispatchObject
#from .._FunctPtr cimport FunctPtr

cdef extern from "Action_Dihedral.h": 
    cdef cppclass _Action_Dihedral "Action_Dihedral" (_Action):
        _Action_Dihedral() 
        _DispatchObject * Alloc() 
        void Help() 


cdef class Action_Dihedral(Action):
    cdef _Action_Dihedral* thisptr

