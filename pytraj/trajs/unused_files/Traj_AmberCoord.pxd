# distutils: language = c++
from pytraj.trajs.TrajectoryIO cimport _TrajectoryIO, TrajectoryIO
from pytraj.BufferedFrame cimport _BufferedFrame, BufferedFrame
from pytraj.BaseIOtype cimport _BaseIOtype, BaseIOtype


cdef extern from "Traj_AmberCoord.h": 
    cdef cppclass _Traj_AmberCoord "Traj_AmberCoord" (_TrajectoryIO):
        _Traj_AmberCoord() 
        _BaseIOtype * Alloc() 
        void WriteHelp()

cdef class Traj_AmberCoord (TrajectoryIO):
    cdef _Traj_AmberCoord* thisptr
