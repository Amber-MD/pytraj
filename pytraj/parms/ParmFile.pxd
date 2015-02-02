# distutils: language = c++
from libcpp.string cimport string
from pytraj.parms.ParmIO cimport _ParmIO, ParmIO
from pytraj.FileTypes cimport _FileTypes, FileTypes
from pytraj.FileName cimport _FileName, FileName
from pytraj.Topology cimport _Topology, Topology
from pytraj.ArgList cimport _ArgList, ArgList


cdef extern from "ParmFile.h": 
    ctypedef enum ParmFormatType "ParmFile::ParmFormatType":
        pass
    # we dont need enum here since we keep them in cpptraj_dict.pyx
    #    AMBERPARM "ParmFile::AMBERPARM"
    #    PDBFILE "ParmFile::PDBFILE"
    #    MOL2FILE "ParmFile::MOL2FILE"
    #    CHARMMPSF "ParmFile::CHARMMPSF"
    #    CIFFILE "ParmFile::CIFFILE"
    #    SDFFILE "ParmFile::SDFFILE"
    # use UNKNOWN_PARM as default
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

