# distutils: language = c++
from pytraj.analyses.Analysis cimport _Analysis, Analysis, RetType
from pytraj.core.DispatchObject cimport _DispatchObject, DispatchObject
from pytraj.core._FunctPtr cimport FunctPtr


cdef extern from "Analysis_Modes.h": 
    cdef cppclass _Analysis_Modes "Analysis_Modes" (_Analysis):
        _Analysis_Modes() 
        _DispatchObject * Alloc() 
        void Help()


cdef class Analysis_Modes (Analysis):
    cdef _Analysis_Modes* thisptr

