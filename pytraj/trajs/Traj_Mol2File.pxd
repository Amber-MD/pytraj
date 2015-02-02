# distutils: language = c++
from libcpp.string cimport string
from TrajectoryIO cimport *
from Mol2File cimport *


cdef extern from "Traj_Mol2File.h": 
    # Traj_Mol2File.h
    ctypedef enum MOL2WRITEMODE "Traj_Mol2File::MOL2WRITEMODE":
        NONE "Traj_Mol2File::NONE"
        SINGLE "Traj_Mol2File::SINGLE"
        MOL "Traj_Mol2File::MOL"
        MULTI "Traj_Mol2File::MULTI"
    cdef cppclass _Traj_Mol2File "Traj_Mol2File":
        _Traj_Mol2File() 
        _BaseIOtype * Alloc() 


cdef class Traj_Mol2File:
    cdef _Traj_Mol2File* thisptr

