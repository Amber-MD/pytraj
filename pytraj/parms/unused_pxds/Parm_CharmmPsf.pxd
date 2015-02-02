# distutils: language = c++
from ParmIO cimport *


cdef extern from "Parm_CharmmPsf.h": 
    cdef cppclass _Parm_CharmmPsf "Parm_CharmmPsf":
        _BaseIOtype * Alloc() 
        bint ID_ParmFormat(_CpptrajFile &)
        int processReadArgs(_ArgList &)
        int ReadParm(const string&, _Topology &)
        int WriteParm(const string&, const _Topology&)
        void SetDebug(int)
        int processWriteArgs(_ArgList &)
