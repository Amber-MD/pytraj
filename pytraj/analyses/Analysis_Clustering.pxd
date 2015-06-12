# distutils: language = c++
from pytraj.analyses.Analysis cimport _Analysis, Analysis, RetType
from pytraj.core.DispatchObject cimport _DispatchObject, DispatchObject
from pytraj.core._FunctPtr cimport FunctPtr


cdef extern from "Analysis_Clustering.h": 
    cdef cppclass _Analysis_Clustering "Analysis_Clustering" (_Analysis):
        _Analysis_Clustering() 
        _DispatchObject * Alloc() 
        void Help()


cdef class Analysis_Clustering (Analysis):
    cdef _Analysis_Clustering* thisptr

