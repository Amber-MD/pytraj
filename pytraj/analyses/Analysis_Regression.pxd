# distutils: language = c++
from pytraj.analyses.Analysis cimport _Analysis, Analysis, RetType
from pytraj.core.DispatchObject cimport _DispatchObject, DispatchObject
from pytraj.core._FunctPtr cimport FunctPtr


cdef extern from "Analysis_Regression.h": 
    cdef cppclass _Analysis_Regression "Analysis_Regression" (_Analysis):
        _Analysis_Regression() 
        _DispatchObject * Alloc() 
        void Help()


cdef class Analysis_Regression (Analysis):
    cdef _Analysis_Regression* thisptr

