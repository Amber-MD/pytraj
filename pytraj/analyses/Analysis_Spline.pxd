# distutils: language = c++
from pytraj.analyses.Analysis cimport _Analysis, Analysis, RetType
from pytraj.DispatchObject cimport _DispatchObject, DispatchObject
from pytraj.ArgList cimport _ArgList, ArgList
from pytraj.core.DataFileList cimport _DataFileList, DataFileList
from pytraj.datasets.DataSetList cimport _DataSetList, DataSetList
from pytraj.TopologyList cimport _TopologyList, TopologyList
from pytraj._FunctPtr cimport FunctPtr


cdef extern from "Analysis_Spline.h": 
    cdef cppclass _Analysis_Spline "Analysis_Spline" (_Analysis):
        _Analysis_Spline() 
        _DispatchObject * Alloc() 
        void Help() 
        RetType Setup(_ArgList&, _DataSetList *, _TopologyList *, _DataFileList *, int)
        RetType Analyze() 


cdef class Analysis_Spline (Analysis):
    cdef _Analysis_Spline* thisptr

