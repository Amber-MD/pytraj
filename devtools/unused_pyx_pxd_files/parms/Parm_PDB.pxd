# distutils: language = c++
from ParmIO cimport *


cdef extern from "Parm_PDB.h": 
    cdef cppclass _Parm_PDB "Parm_PDB":
        _Parm_PDB() : readAsPQR_(false)
        _BaseIOtype * Alloc() 
        void ReadHelp() 
        bint ID_ParmFormat(_CpptrajFile &)
        int processReadArgs(_ArgList &)
        int ReadParm(const string&, _Topology &)
        int WriteParm(const string&, const _Topology&)
        void SetDebug(int)
        int processWriteArgs(_ArgList &)
