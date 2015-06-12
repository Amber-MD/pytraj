# distutils: language = c++
from pytraj.analyses.Analysis cimport _Analysis, Analysis, RetType
from pytraj.core.DispatchObject cimport _DispatchObject, DispatchObject
from pytraj.core._FunctPtr cimport FunctPtr


cdef extern from "Analysis_Timecorr.h": 
    cdef cppclass _Analysis_Timecorr "Analysis_Timecorr" (_Analysis):
        _Analysis_Timecorr() 
        _DispatchObject * Alloc() 
        void Help()


cdef class Analysis_Timecorr (Analysis):
    cdef _Analysis_Timecorr* thisptr

