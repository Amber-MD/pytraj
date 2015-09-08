# distutil: language = c++

from libcpp.string cimport string
from .cpptraj_core cimport _DispatchObject, DispatchObject, DispatchAllocatorType, FunctPtr
from ..datafiles.datafiles cimport _DataFileList, DataFileList
from ..Topology cimport _Topology, Topology, _TopologyList, TopologyList
from ..Frame cimport _Frame, Frame
from .cpptraj_core cimport _ArgList, ArgList, _AtomMask, AtomMask
from ..datasets.DataSetList cimport _DataSetList, DataSetList

cdef extern from "ActionList.h":
    cdef cppclass _ActionList "ActionList":
        _ActionList()
        void Clear()
        void SetDebug(int)
        int Debug()
        int AddAction(DispatchAllocatorType, _ArgList&,
                      _TopologyList*, 
                      _DataSetList*, _DataFileList*)
        int SetupActions(_Topology**)
        bint DoActions(_Frame **, int)
        void Print()
        void List()
        bint Empty()
        int Naction()
        const string& CmdString(int)
        DispatchAllocatorType ActionAlloc(int i)

cdef class ActionList:
    cdef _ActionList* thisptr

    # alias for TopologyList (self.process(top))
    cdef object toplist

    # check if self.process is already called or not
    cdef bint top_is_processed
