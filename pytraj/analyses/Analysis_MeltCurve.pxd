# distutils: language = c++
from pytraj.analyses.Analysis cimport _Analysis, Analysis, RetType
from pytraj.core.DispatchObject cimport _DispatchObject, DispatchObject
from pytraj.core._FunctPtr cimport FunctPtr


cdef extern from "Analysis_MeltCurve.h": 
    cdef cppclass _Analysis_MeltCurve "Analysis_MeltCurve" (_Analysis):
        _Analysis_MeltCurve() 
        _DispatchObject * Alloc() 
        void Help()


cdef class Analysis_MeltCurve (Analysis):
    cdef _Analysis_MeltCurve* thisptr

