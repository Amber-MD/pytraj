# distutils: language = c++
from pytraj._utils cimport list_to_char_pp
from libcpp.vector cimport vector
from libcpp.string cimport string

cdef class Cpptraj:
    def __cinit__(self):
        self.thisptr = new _Cpptraj()

    def __dealloc__(self):
        del self.thisptr

    def run(self, vector[string] args):
        cdef char** c_argv = list_to_char_pp(args)
        cdef int argc = len(args)
        self.thisptr.RunCpptraj(argc, c_argv)
