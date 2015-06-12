# distutils: language = c++
from libcpp.vector cimport vector
from libcpp.string cimport string

from .DispatchObject cimport DispatchAllocatorType
from .TopologyList cimport _TopologyList, TopologyList
from .DataFileList cimport _DataFileList, DataFileList
from .ActionList cimport _ActionList, ActionList
from ._FunctPtr cimport FunctPtr
from .DataFile cimport _DataFile, DataFile
from .FileName cimport _FileName, FileName
from .TrajinList cimport _TrajinList, TrajinList

from ..ArgList cimport _ArgList, ArgList
from ..datasets.DataSetList cimport _DataSetList, DataSetList


cdef extern from "CpptrajState.h": 
    cdef cppclass _CpptrajState "CpptrajState":
        _CpptrajState()
        _TopologyList * PFL()
        _DataSetList * DSL()
        _DataFileList * DFL()
        void SetNoExitOnError()
        void SetNoProgress()
        void SetActionSilence(bint b)
        bint ExitOnError()const 
        bint EmptyState()const 
        int AddTrajin(_ArgList &, bint)
        int AddTrajin(const string&)
        int RunAnalyses()
        _TrajinList& InputTrajList()const 
        inline int AddTrajout(const _ArgList&)
        inline int AddTrajout(const string&)
        int AddReference(const string&, _ArgList &)
        inline int AddReference(const string&)
        inline int AddAction(DispatchAllocatorType, _ArgList &)
        inline int AddAnalysis(DispatchAllocatorType, _ArgList &)
        int WorldSize()
        int ListAll(_ArgList &)const 
        int SetListDebug(_ArgList &)
        int ClearList(_ArgList &)
        int RemoveDataSet(_ArgList &)
        int TrajLength(const string&, const vector[string]&)
        int Run()
        void MasterDataFileWrite()

cdef class CpptrajState:
    cdef _CpptrajState* thisptr
    cdef public TopologyList toplist
    cdef public DataFileList datafilelist
    cdef public DataSetList datasetlist
