# distutils: language = c++
from pytraj.analyses.Analysis cimport _Analysis, Analysis, RetType
from pytraj.core.DispatchObject cimport _DispatchObject, DispatchObject
from pytraj.core._FunctPtr cimport FunctPtr


cdef extern from "Analysis_Matrix.h": 
    cdef cppclass _Analysis_Matrix "Analysis_Matrix" (_Analysis):
        _Analysis_Matrix() 
        _DispatchObject * Alloc() 
        void Help()


cdef class Analysis_Matrix (Analysis):
    cdef _Analysis_Matrix* thisptr

