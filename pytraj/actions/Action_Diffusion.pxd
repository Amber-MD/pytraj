# distutils: language = c++
#from libcpp.vector cimport vector
from pytraj.actions.Action cimport _Action, Action, FunctPtr, _DispatchObject


cdef extern from "Action_Diffusion.h": 
    cdef cppclass _Action_Diffusion "Action_Diffusion" (_Action):
        _Action_Diffusion() 
        _DispatchObject * Alloc() 
        void Help() 


cdef class Action_Diffusion (Action):
    cdef _Action_Diffusion* thisptr

