# distutils: language = c++
from libcpp.string cimport string
from ParmIO cimport *


cdef extern from "Parm_Mol2.h": 
    cdef cppclass _Parm_Mol2 "Parm_Mol2":
        _BaseIOtype * Alloc() 
        bint ID_ParmFormat(_CpptrajFile&)
        int processReadArgs(_ArgList&)
        int ReadParm(const string&, _Topology&)
        int WriteParm(const string&, const _Topology&)
        void SetDebug(int)
        int processWriteArgs(_ArgList&)
