# distutils: language = c++
#from libcpp.string cimport string
from pytraj.actions.Action cimport _Action, Action, FunctPtr, _DispatchObject
#from Traj_AmberNetcdf cimport *
#from DataSet_1D cimport *


cdef extern from "Action_CreateReservoir.h": 
    cdef cppclass _Action_CreateReservoir "Action_CreateReservoir" (_Action):
        _Action_CreateReservoir() 
        _DispatchObject * Alloc() 
        void Help() 


cdef class Action_CreateReservoir (Action):
    cdef _Action_CreateReservoir* thisptr

