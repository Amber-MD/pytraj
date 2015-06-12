# distutils: language = c++
from pytraj.analyses.Analysis cimport _Analysis, Analysis, RetType
from pytraj.core.DispatchObject cimport _DispatchObject, DispatchObject
from pytraj.core._FunctPtr cimport FunctPtr


cdef extern from "Analysis_Rms2d.h": 
    cdef cppclass _Analysis_Rms2d "Analysis_Rms2d" (_Analysis):
        _Analysis_Rms2d() 
        _DispatchObject * Alloc() 
        void Help()


cdef class Analysis_Rms2d (Analysis):
    cdef _Analysis_Rms2d* thisptr

