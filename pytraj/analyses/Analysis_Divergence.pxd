# distutils: language = c++
from pytraj.analyses.Analysis cimport _Analysis, Analysis, RetType
from pytraj.core.DispatchObject cimport _DispatchObject, DispatchObject
from pytraj.core._FunctPtr cimport FunctPtr


cdef extern from "Analysis_Divergence.h": 
    cdef cppclass _Analysis_Divergence "Analysis_Divergence" (_Analysis):
        _Analysis_Divergence() 
        _DispatchObject * Alloc() 
        void Help()


cdef class Analysis_Divergence (Analysis):
    cdef _Analysis_Divergence* thisptr

