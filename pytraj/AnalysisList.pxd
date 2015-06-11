# distutils: language = c++
from .DispatchObject  cimport DispatchAllocatorType
from .analyses.Analysis cimport _Analysis, Analysis
from .ArgList cimport _ArgList 
from .TopologyList cimport _TopologyList
from .datasets.DataSetList cimport _DataSetList
from .core.DataFileList cimport _DataFileList


cdef extern from "AnalysisList.h": 
    cdef cppclass _AnalysisList "AnalysisList":
        _AnalysisList() 
        #~_AnalysisList() 
        void Clear() 
        void SetDebug(int)
        int Debug() const 
        int AddAnalysis(DispatchAllocatorType, _ArgList &, _TopologyList *, _DataSetList *, _DataFileList *)
        int DoAnalyses() 
        void List() const 
        bint Empty() const 

cdef class AnalysisList:
    cdef _AnalysisList* thisptr
