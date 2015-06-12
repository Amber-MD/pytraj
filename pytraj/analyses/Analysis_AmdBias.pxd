# distutils: language = c++
from pytraj.analyses.Analysis cimport _Analysis, Analysis, RetType
from pytraj.core.DispatchObject cimport _DispatchObject, DispatchObject
from pytraj.core._FunctPtr cimport FunctPtr


cdef extern from "Analysis_AmdBias.h": 
    cdef cppclass _Analysis_AmdBias "Analysis_AmdBias" (_Analysis):
        _Analysis_AmdBias() 
        _DispatchObject * Alloc() 
        void Help()


cdef class Analysis_AmdBias (Analysis):
    cdef _Analysis_AmdBias* thisptr

