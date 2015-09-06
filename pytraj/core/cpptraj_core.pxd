# distutil: language = c++
from libcpp.string cimport string
from pytraj.core.Box cimport _Box, Box
from posix.unistd cimport off_t
from libcpp.string cimport string

ctypedef _BaseIOtype* (*AllocatorType)()
ctypedef void (*HelpType)()

cdef extern from "BaseIOtype.h":
    #ctypedef _BaseIOtype* (*AllocatorType)()
    #ctypedef void (*HelpType)()
    cdef cppclass _BaseIOtype "BaseIOtype":
        pass

cdef class BaseIOtype:
    cdef _BaseIOtype* baseptr0

ctypedef _DispatchObject* (*DispatchAllocatorType)()
cdef extern from "DispatchObject.h":
    cdef cppclass _DispatchObject "DispatchObject":
        pass

cdef class DispatchObject:
    cdef _DispatchObject* thisptr

# dummy class to hold function pointer
cdef class FunctPtr:
    cdef DispatchAllocatorType ptr
    # used for BaseIOtype
    cdef AllocatorType allocptr
# distutils: language = c++

cdef extern from "FileName.h":
    cdef cppclass _FileName "FileName":
        _FileName()
        _FileName(_FileName)
        int SetFileName(string)
        int SetFileNameWithExpansion(string)
        int SetFileName(string, bool)
        void clear()
        bint MatchFullOrBase(string)
        string Full()
        string Base()
        char * full()
        char * base()
        string Ext()
        string Compress()
        string DirPrefix()
        bint empty()

cdef class FileName:
    cdef _FileName* thisptr

cdef extern from "CoordinateInfo.h": 
    cdef cppclass _CoordinateInfo "CoordinateInfo":
        _CoordinateInfo() 
        _CoordinateInfo(const _Box& b, bint v, bint t, bint m)
        bint HasBox() const 
        const _Box& TrajBox() const 
        bint HasVel() const 
        bint HasTemp() const 
        bint HasTime() const 
        bint HasForce() const 
        bint HasReplicaDims() const 
        void SetTime(bint m)
        void SetTemperature(bint t)
        void SetVelocity(bint v)
        void SetBox(const _Box& b)

cdef class CoordinateInfo:
    cdef _CoordinateInfo* thisptr

# distutil: language = c++

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
# distutil: language = c++

cdef extern from "NameType.h":
    cdef cppclass _NameType "NameType":
        _NameType() 
        _NameType(const _NameType&)
        _NameType(const char *)
        _NameType(const string&)
        #_NameType& operator =(const _NameType&)
        void ToBuffer(char *) const 
        bint Match(const _NameType&) const 
        bint operator ==(const _NameType&) const 
        bint operator ==(const char *) const 
        #bint opr_ne "operator !="(const _NameType&) const 
        bint operator !=(const _NameType&) const 
        bint operator !=(const char *) const 
        const char* opr_star "operator*" () const 
        char opr_idx "operator[]"(int) const 
        string Truncated() const 
        void ReplaceAsterisk() 


cdef class NameType:
        cdef _NameType* thisptr
