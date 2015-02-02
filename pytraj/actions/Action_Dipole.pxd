# distutils: language = c++
#from libcpp.vector cimport vector
#from libcpp.string cimport string
from pytraj.actions.Action cimport _Action, Action, FunctPtr, _DispatchObject
#from DataSet_GridFlt cimport *
#from GridAction cimport *


cdef extern from "Action_Dipole.h": 
    cdef cppclass _Action_Dipole "Action_Dipole" (_Action):
        _Action_Dipole() 
        _DispatchObject * Alloc() 
        void Help() 


cdef class Action_Dipole (Action):
    cdef _Action_Dipole* thisptr

