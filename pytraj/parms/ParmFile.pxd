# distutils: language = c++
from libcpp.string cimport string
from ..core.cpptraj_core cimport _FileName, FileName
from ..Topology cimport _Topology, Topology
from ..ArgList cimport _ArgList, ArgList


cdef extern from "ParmFile.h": 
    ctypedef enum ParmFormatType "ParmFile::ParmFormatType":
        pass
        UNKNOWN_PARM "ParmFile::UNKNOWN_PARM"
    cdef cppclass _ParmFile "ParmFile":
        @staticmethod
        void ReadOptions() 
        @staticmethod
        void WriteOptions() 
        _ParmFile() 
        int ReadTopology(_Topology&, const string&, const _ArgList&, int)
        int ReadTopology(_Topology& t, const string& n, int d)
        int WritePrefixTopology(const _Topology&, const string&, ParmFormatType, int)
        int WriteTopology(const _Topology&, const string&, const _ArgList&, ParmFormatType, int)
        int WriteTopology(const _Topology& t, const string& n, ParmFormatType f, int d)
        const _FileName ParmFilename() 


cdef class ParmFile:
    cdef _ParmFile* thisptr
