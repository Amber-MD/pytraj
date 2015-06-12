# distutils: language = c++
from libcpp.string cimport string
from pytraj.ArgList cimport _ArgList, ArgList
from pytraj.Topology cimport Topology, _Topology
from pytraj.CpptrajFile cimport _CpptrajFile, CpptrajFile
from pytraj.BaseIOtype cimport _BaseIOtype, BaseIOtype


cdef extern from "ParmIO.h": 
    cdef cppclass _ParmIO "ParmIO":
        # this is Abstract Base Class
        # below are virtual methods
        #~_ParmIO() 
        bint ID_ParmFormat(_CpptrajFile&) 
        int processReadArgs(_ArgList&) 
        int ReadParm(const string&, _Topology&) 
        int processWriteArgs(_ArgList&) 
        int WriteParm(const string&, const _Topology&) 
        void SetDebug(int) 


cdef class ParmIO:
    # use baseptr (instead of thisptr) as a name for base class 
    cdef _ParmIO* baseptr
