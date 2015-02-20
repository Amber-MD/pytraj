# distutils: language = c++
from pytraj.actions.Action cimport _Action, Action, FunctPtr, _DispatchObject


cdef extern from "Action_AreaPerMol.h": 
    cdef cppclass _Action_AreaPerMol "Action_AreaPerMol" (_Action):
        _Action_AreaPerMol() 
        _DispatchObject * Alloc() 
        void Help() 


cdef class Action_AreaPerMol (Action):
    cdef _Action_AreaPerMol* thisptr

