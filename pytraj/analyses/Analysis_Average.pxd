# distutils: language = c++
from pytraj.analyses.Analysis cimport _Analysis, Analysis
from pytraj.DispatchObject cimport _DispatchObject, DispatchObject
from pytraj.ArgList cimport _ArgList, ArgList
from pytraj.DataFileList cimport _DataFileList, DataFileList
from pytraj.DataSetList cimport _DataSetList, DataSetList
from pytraj.TopologyList cimport _TopologyList, TopologyList
from pytraj._FunctPtr cimport FunctPtr


cdef extern from "Analysis_Average.h": 
    cdef cppclass _Analysis_Average "Analysis_Average" (_Analysis):
        _Analysis_Average() 
        _DispatchObject * Alloc() 
        void Help() 
        RetType Setup(_ArgList&, _DataSetList *, _TopologyList *, _DataFileList *, int)
        RetType Analyze() 


cdef class Analysis_Average (Analysis):
    cdef _Analysis_Average* thisptr

