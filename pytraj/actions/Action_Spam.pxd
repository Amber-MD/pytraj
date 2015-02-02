# distutils: language = c++
#from libcpp.vector cimport vector
#from libcpp.string cimport string
from pytraj.actions.Action cimport _Action, Action, FunctPtr, _DispatchObject
#from ImagedAction cimport *
#from Vec3 cimport *


cdef extern from "Action_Spam.h": 
    cdef cppclass _Action_Spam "Action_Spam" (_Action):
        _Action_Spam() 
        _DispatchObject * Alloc() 
        void Help() 


cdef class Action_Spam (Action):
    cdef _Action_Spam* thisptr

