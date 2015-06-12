# distutils: language = c++
from pytraj.analyses.Analysis cimport _Analysis, Analysis, RetType
from pytraj.core.DispatchObject cimport _DispatchObject, DispatchObject
from pytraj.core._FunctPtr cimport FunctPtr


cdef extern from "Analysis_FFT.h": 
    cdef cppclass _Analysis_FFT "Analysis_FFT" (_Analysis):
        _Analysis_FFT() 
        _DispatchObject * Alloc() 
        void Help()


cdef class Analysis_FFT (Analysis):
    cdef _Analysis_FFT* thisptr

