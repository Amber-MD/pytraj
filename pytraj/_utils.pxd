# distutils: language = c++
from libc.stdlib cimport malloc 
from libcpp.vector cimport vector
from libcpp.string cimport string
from cython cimport view
cimport cython
from cython.view cimport array as cyarray
from cpython.array import array as pyarray

num_fused_type = cython.fused_type(int, float, double)

cdef inline int get_positive_idx(idx, size):
    # TODO : do we need this method?
    # we can we memoryview to get slicing too
    """Used for negative indexing"""
    if idx < 0:
        idx = size + idx
        if idx < 0:
            raise ValueError("index is out of range")
    if idx >= size:
        raise ValueError("index is out of range")
    return idx
