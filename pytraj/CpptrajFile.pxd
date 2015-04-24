# distutil: language = c++

from posix.unistd cimport off_t
from libcpp.string cimport string
from .FileName cimport *
#from pytraj.FileIO  cimport *

cdef extern from "CpptrajFile.h":
    ctypedef enum AccessType "CpptrajFile::AccessType":
        pass
    ctypedef enum CompressType "CpptrajFile::CompressType":
        pass
    ctypedef enum FileType "CpptrajFile::FileType":
        pass
    cdef cppclass _CpptrajFile "CpptrajFile":
        _CpptrajFile()
        _CpptrajFile(const _CpptrajFile&)
        int OpenRead(const string&)
        int SetupRead(const string&, int)
        int OpenWriteNumbered(int)
        int OpenWrite(const string&)
        #int OpenEnsembleWrite(const string&, int)
        int SetupWrite(const string&, int)
        int SetupWrite(const string&, FileType, int)
        int OpenAppend(const string&)
        #int OpenEnsembleAppend(const string&, int)
        int SetupAppend(const string&, int)
        int OpenFile()
        int OpenFile(AccessType)
        void CloseFile()
        void Printf(const char*, ...)
        string GetLine() except +
        const char* NextLine()
        AccessType Access()
        CompressType Compression()
        bint IsOpen()
        const _FileName& Filename()
        int IsDos()
        off_t FileSize()
        bint IsCompressed()
        off_t UncompressedSize()
        int ID_Type(const char* filenameIn)
        #int Gets(char*, int)
        #int Write(const void*, size_t)
        int Read(void*, size_t)
        #int Seek(off_t)
        #int Rewind()
        #int Flush()
        #off_t Tell()

cdef class CpptrajFile:
    cdef _CpptrajFile* thisptr
