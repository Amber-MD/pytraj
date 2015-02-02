# distutils: language = c++
#from libcpp.vector cimport vector
from pytraj.actions.Action cimport _Action, Action, FunctPtr, _DispatchObject
#from Trajout cimport *
#from DataSet_1D cimport *


cdef extern from "Action_Outtraj.h": 
    cdef cppclass _Action_Outtraj "Action_Outtraj" (_Action):
        _Action_Outtraj() 
        _DispatchObject * Alloc() 
        void Help() 


cdef class Action_Outtraj (Action):
    cdef _Action_Outtraj* thisptr

