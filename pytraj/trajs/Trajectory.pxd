# distutils: language = c++
from libcpp.string cimport string
from libcpp.vector cimport vector
from ..Frame cimport _Frame, Frame
from ..Topology cimport _Topology, Topology
#from .trajs.Trajin_Single cimport _Trajin_Single, Trajin_Single
from .Trajin cimport _Trajin, Trajin
from ..cpp_algorithm cimport reverse as cpp_reverse

ctypedef vector[_Frame*].iterator iterator


cdef class Trajectory:
    cdef vector[_Frame*] frame_v
    cdef public Topology top
    cdef public Topology oldtop 

    # used for warning memory view
    cdef public bint warning 
    #cdef void _join(Trajectory self, Trajectory other, mask=*)

    # create tmpfarray to hold sub Trajectory
    # traj[0:10][0] will give wrong answer
    cdef object tmpfarray

    # make public?
    # this variable is intended to let Trajectory control 
    # freeing memory for Frame instance but it's too complicated
    # 
    #cdef bint is_mem_parent
