# distutils: language = c++
#from libcpp.vector cimport vector
#from libcpp.string cimport string
from pytraj.actions.Action cimport _Action, Action, FunctPtr, _DispatchObject
#from DihedralSearch cimport *


cdef extern from "Action_MakeStructure.h": 
#    cdef cppclass _Action_MakeStructure::SS_TYPE "Action_MakeStructure::SS_TYPE" (_Action):
#        SS_TYPE() 
#        SS_TYPE(double ph, double ps, double ph2, double ps2, int t, const string& n)
#        bint empty() 
#
#
#    cdef cppclass _Action_MakeStructure::SecStructHolder "Action_MakeStructure::SecStructHolder" (_Action):
#        SecStructHolder() 
#        SecStructHolder(const string& rangearg, int typeidx)


    cdef cppclass _Action_MakeStructure "Action_MakeStructure" (_Action):
        _Action_MakeStructure() 
        _DispatchObject * Alloc() 
        void Help() 

cdef class Action_MakeStructure (Action):
    cdef _Action_MakeStructure* thisptr

