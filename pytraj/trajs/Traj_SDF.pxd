# distutils: language = c++
from TrajectoryIO cimport *
from SDFfile cimport *


cdef extern from "Traj_SDF.h": 
    cdef cppclass _Traj_SDF "Traj_SDF":
        _Traj_SDF() 
        _BaseIOtype * Alloc() 
