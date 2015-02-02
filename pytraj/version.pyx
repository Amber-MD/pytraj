# distutils: language = c++
from libcpp.string cimport string

cdef extern from "Version.h":
    cdef string CPPTRAJ_VERSION_STRING

cppversion = CPPTRAJ_VERSION_STRING.decode()
