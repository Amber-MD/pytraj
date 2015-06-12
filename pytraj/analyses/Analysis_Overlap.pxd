# distutils: language = c++
from pytraj.analyses.Analysis cimport _Analysis, Analysis, RetType
from pytraj.core.DispatchObject cimport _DispatchObject, DispatchObject
from pytraj.core._FunctPtr cimport FunctPtr


cdef extern from "Analysis_Overlap.h": 
    cdef cppclass _Analysis_Overlap "Analysis_Overlap" (_Analysis):
        _Analysis_Overlap() 
        _DispatchObject * Alloc() 
        void Help()


cdef class Analysis_Overlap (Analysis):
    cdef _Analysis_Overlap* thisptr

