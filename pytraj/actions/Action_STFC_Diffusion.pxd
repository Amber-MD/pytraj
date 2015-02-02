# distutils: language = c++
#from libcpp.vector cimport vector
#from libcpp.string cimport string
from pytraj.actions.Action cimport _Action, Action, FunctPtr, _DispatchObject


cdef extern from "Action_STFC_Diffusion.h": 
    cdef cppclass _Action_STFC_Diffusion "Action_STFC_Diffusion" (_Action):
        _Action_STFC_Diffusion() 
        _DispatchObject * Alloc() 
        void Help() 


cdef class Action_STFC_Diffusion (Action):
    cdef _Action_STFC_Diffusion* thisptr

