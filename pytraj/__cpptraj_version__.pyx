# distutils: language = c++
from libcpp.string cimport string
from libcpp.vector cimport vector

cdef extern from "Version.h":
    cdef string CPPTRAJ_VERSION_STRING
    cdef string CPPTRAJ_INTERNAL_VERSION

cdef extern from "Cpptraj.h":
    cdef cppclass Cpptraj:
        @staticmethod
        string Defines()

__cpptraj_version__ = CPPTRAJ_VERSION_STRING.decode()
__cpptraj_internal_version__ = CPPTRAJ_INTERNAL_VERSION.decode()

def info():
    cdef string s

    s = Cpptraj.Defines()
    return s.decode()
