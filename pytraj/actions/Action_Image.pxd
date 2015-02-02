# distutils: language = c++
from pytraj.actions.Action cimport _Action, Action, FunctPtr, _DispatchObject
#from ImageTypes cimport *


cdef extern from "Action_Image.h": 
    cdef cppclass _Action_Image "Action_Image" (_Action):
        _Action_Image() 
        _DispatchObject * Alloc() 
        void Help() 
        #~_Action_Image() 


cdef class Action_Image (Action):
    cdef _Action_Image* thisptr

