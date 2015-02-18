# distutils: language = c++
from TrajectoryIO cimport *
from CIFfile cimport *


cdef extern from "Traj_CIF.h": 
    cdef cppclass _Traj_CIF "Traj_CIF":
        _Traj_CIF() : Natoms_(0 ), Nmodels_(0 ), Cartn_x_col_(0 ), Cartn_y_col_(0 ), Cartn_z_col_(0)
        _BaseIOtype * Alloc() 
