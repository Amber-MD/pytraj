# distutils: language = c++
from pytraj.analyses.Analysis cimport _Analysis, Analysis, RetType
from pytraj.core.DispatchObject cimport _DispatchObject, DispatchObject
from pytraj.core._FunctPtr cimport FunctPtr


cdef extern from "Analysis_IRED.h": 
    cdef cppclass _Analysis_IRED "Analysis_IRED" (_Analysis):
        _Analysis_IRED() 
        _DispatchObject * Alloc() 
        void Help()


cdef class Analysis_IRED (Analysis):
    cdef _Analysis_IRED* thisptr

