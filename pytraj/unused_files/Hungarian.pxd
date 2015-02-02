# distutils: language = c++
from libcpp.vector cimport vector
from pycpptraj.Matrix cimport *


cdef extern from "Hungarian.h": 
    cdef cppclass _Hungarian "Hungarian":
        _Hungarian()
        int Initialize(size_t)
        void AddElement(double d)
        vector [int] Optimize() 
