# distutils: language = c++
from __future__ import absolute_import
from libcpp.string cimport string
from ..ArgList cimport _ArgList, ArgList
from ..core.FileName cimport FileName, _FileName
from ..Topology cimport _Topology, Topology


cdef extern from "TrajectoryFile.h": 
    ctypedef enum TrajFormatType "TrajectoryFile::TrajFormatType":
        pass

    cdef cppclass _TrajectoryFile "TrajectoryFile":
        _TrajectoryFile ()
        @staticmethod
        void ReadOptions ()
        @staticmethod
        void WriteOptions ()
        @staticmethod
        TrajFormatType GetFormatFromArg(_ArgList & a)
        @staticmethod
        TrajFormatType GetFormatFromString(const string& s)
        @staticmethod
        string GetExtensionForType(TrajFormatType t)
        @staticmethod
        TrajFormatType GetTypeFromExtension(const string& e)
        @staticmethod
        const char * FormatString(TrajFormatType tt)
        void SetDebug(int)
        void SetTrajFileName(const string&, bint)
        int SetTrajParm(_Topology *)
        _Topology * TrajParm ()const 
        const _FileName & TrajFilename ()const 


cdef class TrajectoryFile:
    cdef _TrajectoryFile* baseptr0
    cdef Topology _top
