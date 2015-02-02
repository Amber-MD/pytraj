# distutils: language = c++
from libcpp.string cimport string
from pytraj.parms.ParmIO cimport _ParmIO, ParmIO
from pytraj.ArgList cimport _ArgList, ArgList
from pytraj.CpptrajFile cimport _CpptrajFile, CpptrajFile
from pytraj.BaseIOtype cimport _BaseIOtype, BaseIOtype
from pytraj.Topology cimport _Topology, Topology


cdef extern from "Parm_Amber.h": 
    cdef cppclass _Parm_Amber "Parm_Amber":
        _Parm_Amber() 
        #~_Parm_Amber() 
        _BaseIOtype * Alloc() 
        void WriteHelp() 
        bint ID_ParmFormat(_CpptrajFile&)
        int processReadArgs(_ArgList&)
        int ReadParm(const string&, _Topology&)
        int processWriteArgs(_ArgList&)
        int WriteParm(const string&, const _Topology&)
        void SetDebug(int debugIn)

cdef class Parm_Amber:
    cdef _Parm_Amber* thisptr
