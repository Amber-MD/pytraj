# distutils: language = c++
from libcpp.string cimport string
from libcpp.vector cimport vector
from ..Frame cimport _Frame, Frame
from ..Topology cimport _Topology, Topology
from ..cpp_algorithm cimport reverse as cpp_reverse

ctypedef vector[_Frame*].iterator iterator


cdef class Trajectory:
    cdef vector[_Frame*] frame_v
    cdef public Topology top
    cdef public Topology oldtop 
    # used for warning memory view
    cdef public bint warning 
    # create tmpfarray to hold sub Trajectory
    cdef object tmpfarray
