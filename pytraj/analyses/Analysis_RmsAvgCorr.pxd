# distutils: language = c++
from pytraj.analyses.Analysis cimport _Analysis, Analysis, RetType
from pytraj.core.DispatchObject cimport _DispatchObject, DispatchObject
from pytraj.core._FunctPtr cimport FunctPtr


cdef extern from "Analysis_RmsAvgCorr.h": 
    cdef cppclass _Analysis_RmsAvgCorr "Analysis_RmsAvgCorr" (_Analysis):
        _Analysis_RmsAvgCorr() 
        _DispatchObject * Alloc() 
        void Help()


cdef class Analysis_RmsAvgCorr (Analysis):
    cdef _Analysis_RmsAvgCorr* thisptr

