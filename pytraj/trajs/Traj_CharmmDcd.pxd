# distutils: language = c++
from TrajectoryIO cimport *


cdef extern from "Traj_CharmmDcd.h": 
    cdef cppclass _Traj_CharmmDcd "Traj_CharmmDcd":
        _Traj_CharmmDcd() 
        _BaseIOtype * Alloc() 
        #~_Traj_CharmmDcd() 
