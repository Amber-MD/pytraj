# distutils: language = c++
from pytraj.DispatchObject cimport DispatchAllocatorType 
from pytraj.BaseIOtype cimport AllocatorType
from pytraj.trajs.Trajin cimport _Trajin

# dummy class to hold function pointer
cdef class FunctPtr:
    cdef DispatchAllocatorType ptr
    # used for BaseIOtype
    cdef AllocatorType allocptr
    cdef _Trajin* trajinptr
