# distutil: language = c++

from libcpp.string cimport string
from pytraj.DispatchObject cimport _DispatchObject, DispatchObject, DispatchAllocatorType
from pytraj.Topology cimport _Topology, Topology
from pytraj.TopologyList cimport _TopologyList, TopologyList
from pytraj.Frame cimport _Frame, Frame
from pytraj.DataSetList cimport _DataSetList, DataSetList
from pytraj.DataFileList cimport _DataFileList, DataFileList
from pytraj.FrameArray cimport FrameArray
from pytraj.ArgList cimport _ArgList, ArgList
from pytraj._FunctPtr cimport FunctPtr
from pytraj.AtomMask cimport _AtomMask, AtomMask
#from pytraj.CpptrajFile cimport _CpptrajFile, CpptrajFile

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
