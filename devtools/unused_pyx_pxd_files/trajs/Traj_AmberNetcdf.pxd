# distutils: language = c++
from libcpp.string cimport string
from TrajectoryIO cimport *
from NetcdfFile cimport *
include "config.pxi"


cdef extern from "Traj_AmberNetcdf.h": 
    cdef cppclass _Traj_AmberNetcdf "Traj_AmberNetcdf":
        _Traj_AmberNetcdf() 
        _BaseIOtype * Alloc() 
        #~_Traj_AmberNetcdf() 
        bint ID_TrajFormat(_CpptrajFile&)
        int setupTrajin(const string&, _Topology *)
        int setupTrajout(const string&, _Topology *, int, bint)
        int openTrajin() 
        void closeTraj() 
        int readFrame(int, _Frame&)
        int readVelocity(int, _Frame&)
        int writeFrame(int, const _Frame&)
        void Info() 
        int processWriteArgs(_ArgList&)
        int processReadArgs(_ArgList&)
        inline int createReservoir(bint, double, int)
        int writeReservoir(int, _Frame&, double, int)


IF BINTRAJ: 
    cdef class Traj_AmberNetcdf:
        cdef _Traj_AmberNetcdf* thisptr

