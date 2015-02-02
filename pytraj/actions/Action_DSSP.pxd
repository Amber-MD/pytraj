# distutils: language = c++
#from libcpp.vector cimport vector
#from libcpp.string cimport string
from pytraj.actions.Action cimport _Action, Action, FunctPtr, _DispatchObject


cdef extern from "Action_DSSP.h": 
    cdef cppclass _Action_DSSP "Action_DSSP" (_Action):
        _Action_DSSP() 
        _DispatchObject * Alloc() 
        void Help() 


cdef class Action_DSSP (Action):
    cdef _Action_DSSP* thisptr

