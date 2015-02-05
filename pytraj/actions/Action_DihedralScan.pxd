# distutils: language = c++
#from libcpp.vector cimport vector
#from libcpp.string cimport string
from pytraj.actions.Action cimport _Action, Action, FunctPtr, _DispatchObject
#from Action_CheckStructure cimport *
#from Random cimport *
#from Trajout cimport *
#from DihedralSearch cimport *


cdef extern from "Action_DihedralScan.h": 
    cdef cppclass _Action_DihedralScan "Action_DihedralScan" (_Action):
        _Action_DihedralScan() 
        #~_Action_DihedralScan() 
        _DispatchObject * Alloc() 
        void Help() 

cdef class Action_DihedralScan (Action):
    cdef _Action_DihedralScan* thisptr

