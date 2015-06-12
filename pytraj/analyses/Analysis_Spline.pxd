# distutils: language = c++
from pytraj.analyses.Analysis cimport _Analysis, Analysis, RetType
from pytraj.core.DispatchObject cimport _DispatchObject, DispatchObject
from pytraj.core._FunctPtr cimport FunctPtr


cdef extern from "Analysis_Spline.h": 
    cdef cppclass _Analysis_Spline "Analysis_Spline" (_Analysis):
        _Analysis_Spline() 
        _DispatchObject * Alloc() 
        void Help()


cdef class Analysis_Spline (Analysis):
    cdef _Analysis_Spline* thisptr

