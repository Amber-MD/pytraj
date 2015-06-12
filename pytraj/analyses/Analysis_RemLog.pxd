# distutils: language = c++
from pytraj.analyses.Analysis cimport _Analysis, Analysis, RetType
from pytraj.core.DispatchObject cimport _DispatchObject, DispatchObject
from pytraj.core._FunctPtr cimport FunctPtr


cdef extern from "Analysis_RemLog.h": 
    cdef cppclass _Analysis_RemLog "Analysis_RemLog" (_Analysis):
        _Analysis_RemLog() 
        _DispatchObject * Alloc() 
        void Help()


cdef class Analysis_RemLog (Analysis):
    cdef _Analysis_RemLog* thisptr

