# distutils: language = c++
#from libcpp.vector cimport vector
from pytraj.actions.Action cimport _Action, Action, FunctPtr, _DispatchObject


cdef extern from "Action_Contacts.h": 
    cdef cppclass _Action_Contacts "Action_Contacts" (_Action):
        _Action_Contacts() 
        _DispatchObject * Alloc() 
        void Help() 
        #~_Action_Contacts() 


cdef class Action_Contacts (Action):
    cdef _Action_Contacts* thisptr

