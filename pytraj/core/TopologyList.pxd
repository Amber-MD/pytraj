# distutil: language = c++
from libcpp.string cimport string
from pytraj.Topology cimport _Topology, Topology
from pytraj.ArgList cimport _ArgList, ArgList

cdef extern from "TopologyList.h":
    cdef cppclass _TopologyList "TopologyList":
        const char* ParmArgs
        TopologyList()
        void Clear()
        void SetDebug(int)
        _Topology* GetParm(int) 
        _Topology* GetParmByIndex(_ArgList&) 
        _Topology* GetParm(_ArgList&) 
        int AddParmFile(string&)
        int AddParmFile(string&, _ArgList&)
        void AddParm(_Topology * pIn)
        void List()

cdef class TopologyList:
    cdef _TopologyList* thisptr
    cdef bint py_free_mem
