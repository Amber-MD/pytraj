# distutils: language = c++
from pytraj.analyses.Analysis cimport _Analysis, Analysis, RetType
from pytraj.core.DispatchObject cimport _DispatchObject, DispatchObject
from pytraj.core._FunctPtr cimport FunctPtr


cdef extern from "Analysis_Average.h": 
    cdef cppclass _Analysis_Average "Analysis_Average" (_Analysis):
        _Analysis_Average() 
        _DispatchObject * Alloc() 
        void Help()


cdef class Analysis_Average (Analysis):
    cdef _Analysis_Average* thisptr

