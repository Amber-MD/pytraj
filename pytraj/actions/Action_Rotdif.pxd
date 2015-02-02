# distutils: language = c++
#from libcpp.vector cimport vector
#from libcpp.string cimport string
from pytraj.actions.Action cimport _Action, Action, FunctPtr, _DispatchObject
#from Random cimport *
#from DataSet_Vector cimport *


cdef extern from "Action_Rotdif.h": 
    cdef cppclass _Action_Rotdif "Action_Rotdif" (_Action):
        _Action_Rotdif() 
        _DispatchObject * Alloc() 
        void Help() 


cdef class Action_Rotdif (Action):
    cdef _Action_Rotdif* thisptr

