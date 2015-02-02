# this file has several c++ methods
# usage: from pycpptraj.cpp_utils cimport your_desired_code
from libcpp.vector cimport vector

cdef extern from "<std>" namespace "<std>" nogil:
    vector[T] reverse(vector[T].begin(), vector[T].end())
