# distutils: language = c++
from pytraj.DispatchObject  cimport DispatchAllocatorType
from pytraj.analyses.Analysis cimport _Analysis, Analysis
from pytraj.ArgList cimport _ArgList 
from pytraj.TopologyList cimport _TopologyList
from pytraj.DataSetList cimport _DataSetList
from pytraj.DataFileList cimport _DataFileList


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
