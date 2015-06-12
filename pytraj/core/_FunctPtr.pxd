# distutils: language = c++
from .DispatchObject cimport DispatchAllocatorType 
from .BaseIOtype cimport AllocatorType

# dummy class to hold function pointer
cdef class FunctPtr:
    cdef DispatchAllocatorType ptr
    # used for BaseIOtype
    cdef AllocatorType allocptr
