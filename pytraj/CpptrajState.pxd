# distutils: language = c++
from libcpp.vector cimport vector
from libcpp.string cimport string
from pytraj.DispatchObject cimport DispatchAllocatorType
from libcpp.vector cimport vector
from pytraj.TrajinList cimport _TrajinList, TrajinList
from pytraj.ArgList cimport _ArgList, ArgList
from pytraj.TopologyList cimport _TopologyList, TopologyList
from pytraj.DataSetList cimport _DataSetList, DataSetList
from pytraj.DataFileList cimport _DataFileList, DataFileList
from pytraj.ActionList cimport _ActionList, ActionList
from pytraj.AnalysisList cimport _AnalysisList, AnalysisList
from pytraj._FunctPtr cimport FunctPtr
from pytraj.DataFile cimport _DataFile, DataFile
from pytraj.FileName cimport _FileName, FileName


cdef extern from "CpptrajState.h": 
    cdef cppclass _CpptrajState "CpptrajState":
        _CpptrajState()
        _TopologyList * PFL()
        #_FrameList * FL()
        _DataSetList * DSL()
        _DataFileList * DFL()
        void SetNoExitOnError()
        void SetNoProgress()
        void SetActionSilence(bint b)
        #int Debug()const 
        bint ExitOnError()const 
        bint EmptyState()const 
        int AddTrajin(_ArgList &, bint)
        int AddTrajin(const string&)
        int RunAnalyses()
        #const _TrajinList& InputTrajList()const 
        _TrajinList& InputTrajList()const 
        inline int AddTrajout(const _ArgList&)
        inline int AddTrajout(const string&)
        #inline int AddReference(_ArgList &)
        #inline int AddReference(const string&)
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
    #cdef public FrameList framelist
