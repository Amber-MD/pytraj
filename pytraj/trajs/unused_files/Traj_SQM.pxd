# distutils: language = c++
from TrajectoryIO cimport *


cdef extern from "Traj_SQM.h": 
    cdef cppclass _Traj_SQM "Traj_SQM":
        _Traj_SQM() : singleWrite_(false ), chargeIsSet_(false ), charge_(0 ), sqmParm_(0)
        _BaseIOtype * Alloc() 
