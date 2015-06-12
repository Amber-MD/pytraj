# distutils: language = c++
from pytraj.analyses.Analysis cimport _Analysis, Analysis, RetType
from pytraj.core.DispatchObject cimport _DispatchObject, DispatchObject
from pytraj.core._FunctPtr cimport FunctPtr


cdef extern from "Analysis_VectorMath.h": 
    cdef cppclass _Analysis_VectorMath "Analysis_VectorMath" (_Analysis):
        _Analysis_VectorMath() 
        _DispatchObject * Alloc() 
        void Help()


cdef class Analysis_VectorMath (Analysis):
    cdef _Analysis_VectorMath* thisptr

