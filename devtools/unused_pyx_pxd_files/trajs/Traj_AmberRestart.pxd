# distutils: language = c++
from TrajectoryIO cimport *
from BufferedFrame cimport *


cdef extern from "Traj_AmberRestart.h": 
    cdef cppclass _Traj_AmberRestart "Traj_AmberRestart":
        _Traj_AmberRestart() 
        _BaseIOtype * Alloc() 
