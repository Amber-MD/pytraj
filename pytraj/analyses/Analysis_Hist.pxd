# distutils: language = c++
from pytraj.analyses.Analysis cimport _Analysis, Analysis, RetType
from pytraj.core.DispatchObject cimport _DispatchObject, DispatchObject
from pytraj.core._FunctPtr cimport FunctPtr


cdef extern from "Analysis_Hist.h": 
    cdef cppclass _Analysis_Hist "Analysis_Hist" (_Analysis):
        _Analysis_Hist() 
        _DispatchObject * Alloc() 
        void Help()


cdef class Analysis_Hist (Analysis):
    cdef _Analysis_Hist* thisptr

