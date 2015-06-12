# distutils: language = c++
from pytraj.analyses.Analysis cimport _Analysis, Analysis, RetType
from pytraj.core.DispatchObject cimport _DispatchObject, DispatchObject
from pytraj.core._FunctPtr cimport FunctPtr


cdef extern from "Analysis_CrossCorr.h": 
    cdef cppclass _Analysis_CrossCorr "Analysis_CrossCorr" (_Analysis):
        _Analysis_CrossCorr() 
        _DispatchObject * Alloc() 
        void Help()


cdef class Analysis_CrossCorr (Analysis):
    cdef _Analysis_CrossCorr* thisptr

