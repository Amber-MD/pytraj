# distutils: language = c++
from ParmIO cimport *


cdef extern from "Parm_SDF.h": 
    cdef cppclass _Parm_SDF "Parm_SDF":
        _BaseIOtype * Alloc() 
        bint ID_ParmFormat(_CpptrajFile &)
        int processReadArgs(_ArgList &)
        int ReadParm(const string&, _Topology &)
        int WriteParm(const string&, const _Topology&)
        void SetDebug(int)
        int processWriteArgs(_ArgList &)
