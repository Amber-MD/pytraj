# distutils: language = c++
from TrajectoryIO cimport *


cdef extern from "Traj_Conflib.h": 
    cdef cppclass _Traj_Conflib "Traj_Conflib":
        _Traj_Conflib() 
        _BaseIOtype * Alloc() 
