# distutils: language = c++
from pytraj.analyses.Analysis cimport _Analysis, Analysis, RetType
from pytraj.core.DispatchObject cimport _DispatchObject, DispatchObject
from pytraj.core._FunctPtr cimport FunctPtr


cdef extern from "Analysis_AutoCorr.h": 
    cdef cppclass _Analysis_AutoCorr "Analysis_AutoCorr" (_Analysis):
        _Analysis_AutoCorr() 
        _DispatchObject * Alloc() 
        void Help()


cdef class Analysis_AutoCorr (Analysis):
    cdef _Analysis_AutoCorr* thisptr

