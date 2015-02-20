# distutils: language = c++
from pytraj.actions.Action cimport _Action, Action, FunctPtr, _DispatchObject


cdef extern from "Action_MinImage.h": 
    cdef cppclass _Action_MinImage "Action_MinImage" (_Action):
        _Action_MinImage() 
        _DispatchObject * Alloc() 
        void Help() 


cdef class Action_MinImage (Action):
    cdef _Action_MinImage* thisptr

