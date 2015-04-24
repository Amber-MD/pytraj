# distutils: language = c++
from libcpp.vector cimport vector
from libcpp.string cimport string
from pytraj.trajs.TrajectoryIO cimport *
from pytraj.PDBfile cimport *


cdef extern from "Traj_PDBfile.h": 
    # Traj_PDBfile.h
    ctypedef enum PDBWRITEMODE "Traj_PDBfile::PDBWRITEMODE":
        pass
        #NONE "Traj_PDBfile::NONE"
        #SINGLE "Traj_PDBfile::SINGLE"
        #MODEL "Traj_PDBfile::MODEL"
        #MULTI "Traj_PDBfile::MULTI"
    cdef cppclass _Traj_PDBfile "Traj_PDBfile":
        _Traj_PDBfile() 
        _BaseIOtype * Alloc() 
        @staticmethod
        void WriteHelp() 


cdef class Traj_PDBfile:
    cdef _Traj_PDBfile* thisptr

