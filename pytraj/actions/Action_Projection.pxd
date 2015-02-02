# distutils: language = c++
#from libcpp.vector cimport vector
#from pytraj.actions.Action cimport *
from pytraj.actions.Action cimport _Action, Action, FunctPtr, _DispatchObject
#from ..DataSet_Modes cimport *
#from .ActionFrameCounter cimport *
#from ..Array1D cimport *


cdef extern from "Action_Projection.h": 
    cdef cppclass _Action_Projection "Action_Projection" (_Action):
        _Action_Projection() 
        _DispatchObject * Alloc() 
        void Help() 


cdef class Action_Projection (Action):
    cdef _Action_Projection* thisptr

