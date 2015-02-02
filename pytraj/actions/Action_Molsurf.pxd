# distutils: language = c++
#from libcpp.vector cimport vector
from pytraj.actions.Action cimport _Action, Action, FunctPtr, _DispatchObject
#from molsurf cimport *


cdef extern from "Action_Molsurf.h": 
    cdef cppclass _Action_Molsurf "Action_Molsurf" (_Action):
        _Action_Molsurf() 
        _DispatchObject * Alloc() 
        void Help() 
        #~_Action_Molsurf() 


cdef class Action_Molsurf (Action):
    cdef _Action_Molsurf* thisptr

