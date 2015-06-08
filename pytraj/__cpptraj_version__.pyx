# distutils: language = c++
from libcpp.string cimport string
from libcpp.vector cimport vector

cdef extern from "Version.h":
    cdef string CPPTRAJ_VERSION_STRING

cdef extern from "Cpptraj.h":
    cdef cppclass Cpptraj:
        @staticmethod
        string Defines()

__cpptraj_version__ = CPPTRAJ_VERSION_STRING.decode()

def info():
    cdef string s

    s = Cpptraj.Defines()
    return s.decode()
