# distutils: language = c++
from pytraj.analyses.Analysis cimport _Analysis, Analysis, RetType
from pytraj.core.DispatchObject cimport _DispatchObject, DispatchObject
from pytraj.core._FunctPtr cimport FunctPtr


cdef extern from "Analysis_Lifetime.h": 
    cdef cppclass _Analysis_Lifetime "Analysis_Lifetime" (_Analysis):
        _Analysis_Lifetime() 
        _DispatchObject * Alloc() 
        void Help()


cdef class Analysis_Lifetime (Analysis):
    cdef _Analysis_Lifetime* thisptr

