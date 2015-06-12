# distutils: language = c++
from pytraj.analyses.Analysis cimport _Analysis, Analysis, RetType
from pytraj.core.DispatchObject cimport _DispatchObject, DispatchObject
from pytraj.core._FunctPtr cimport FunctPtr


cdef extern from "Analysis_Corr.h": 
    cdef cppclass _Analysis_Corr "Analysis_Corr" (_Analysis):
        _Analysis_Corr() 
        _DispatchObject * Alloc() 
        void Help()


cdef class Analysis_Corr (Analysis):
    cdef _Analysis_Corr* thisptr

