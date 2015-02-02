# distutils: language = c++
#from libcpp.vector cimport vector
#from libcpp.string cimport string
from pytraj.actions.Action cimport _Action, Action, FunctPtr, _DispatchObject
#from DataSet_GridFlt cimport *


cdef extern from "Action_Volmap.h": 
    cdef cppclass _Action_Volmap "Action_Volmap" (_Action):
        _Action_Volmap() 
        _DispatchObject * Alloc() 
        void Help() 


cdef class Action_Volmap (Action):
    cdef _Action_Volmap* thisptr

