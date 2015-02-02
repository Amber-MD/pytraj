# distutil: language = c++

ctypedef _DispatchObject* (*DispatchAllocatorType)()
cdef extern from "DispatchObject.h":
    cdef cppclass _DispatchObject "DispatchObject":
        pass

cdef class DispatchObject:
    cdef _DispatchObject* thisptr
