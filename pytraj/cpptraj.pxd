# distutils: language = c++

from libcpp.vector cimport vector
from libcpp.string cimport string

cdef extern from "Cpptraj.h":
    ctypedef enum Mode "Cpptraj::Mode":
        BATCH "Cpptraj::BATCH"
        ERROR "Cpptraj::ERROR"
        QUIT "Cpptraj::QUIT"
        INTERACTIVE "Cpptraj::INTERACTIVE"
        SILENT_EXIT "Cpptraj::SILENT_EXIT"

    cdef cppclass _Cpptraj "Cpptraj":
        _Cpptraj()
        int RunCpptraj(int, char**)

cdef class Cpptraj:
    cdef _Cpptraj* thisptr
