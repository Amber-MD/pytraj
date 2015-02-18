# distutils: language = c++
from libcpp.string cimport string
from TrajectoryIO cimport *
from NetcdfFile cimport *


cdef extern from "Traj_AmberRestartNC.h": 
    cdef cppclass _Traj_AmberRestartNC "Traj_AmberRestartNC":
        _Traj_AmberRestartNC() 
        _BaseIOtype * Alloc() 
        #~_Traj_AmberRestartNC() 
        bint ID_TrajFormat(_CpptrajFile&)
        int setup_Trajin(const string&, _Topology *)
        int setup_Trajout(const string&, _Topology *, int, bint)
        int open_Trajin() 
        void closeTraj() 
        int read_Frame(int, _Frame&)
        int write_Frame(int, const _Frame&)
        int processWriteArgs(_ArgList&)
        int processReadArgs(_ArgList&)
        void Info() 
