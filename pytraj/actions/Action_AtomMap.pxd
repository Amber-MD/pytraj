# distutils: language = c++
#from libcpp.vector cimport vector
from pytraj.actions.Action cimport _Action, Action, FunctPtr, _DispatchObject

cdef extern from "Action_AtomMap.h": 
    cdef cppclass _Action_AtomMap "Action_AtomMap" (_Action):
        _Action_AtomMap() 
        _DispatchObject * Alloc() 
        void Help() 
        #~_Action_AtomMap() 


cdef class Action_AtomMap (Action):
    cdef _Action_AtomMap* thisptr

