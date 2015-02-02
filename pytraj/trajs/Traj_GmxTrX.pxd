# distutils: language = c++
from TrajectoryIO cimport *


cdef extern from "Traj_GmxTrX.h": 
    cdef cppclass _Traj_GmxTrX "Traj_GmxTrX":
        _Traj_GmxTrX() 
        #~_Traj_GmxTrX() 
        _BaseIOtype * Alloc() 
