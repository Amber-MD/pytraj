# distutils: language = c++
from libcpp.vector cimport vector
from libcpp.string cimport string

from .cpptraj_core cimport DispatchAllocatorType, FunctPtr, _FileName, FileName
from .TopologyList cimport _TopologyList, TopologyList
from .DataFileList cimport _DataFileList, DataFileList
from .DataFile cimport _DataFile, DataFile
from .ActionList cimport _ActionList, ActionList
from .DataFile cimport _DataFile, DataFile

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
        inline int AddTrajout "AddOutputTrajectory" (const _ArgList&)
        inline int AddTrajout "AddOutputTrajectory" (const string&)
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
