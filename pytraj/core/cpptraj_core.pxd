# distutil: language = c++
from __future__ import absolute_import

from libcpp.string cimport string
from posix.unistd cimport off_t
from ..cpp_vector cimport vector

from .Box cimport _Box, Box
from .TopologyList cimport _TopologyList, TopologyList
from ..datafiles.datafiles cimport _DataFileList, DataFileList, _DataFile, DataFile
from ..datasets.DataSetList cimport _DataSetList, DataSetList

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


cdef extern from "AtomMask.h": 
    cdef cppclass _AtomMask "AtomMask":
        _AtomMask()
        _AtomMask(const string&)
        _AtomMask(int, int)
        _AtomMask(int)
        _AtomMask(vector[int], int)
        _AtomMask(const _AtomMask &)
        #_AtomMask & operator =(const _AtomMask &)
        const vector [int]& Selected()const 
        vector[int].const_iterator begin()const 
        vector[int].const_iterator end()const 
        int back()const 
        int Nselected()const 
        const int & index_opr "operator[]"(int idx)const 
        const char * MaskString()const 
        const string& MaskExpression()const 
        bint MaskStringSet()const 
        bint None()const 
        bint IsCharMask()const 
        void ResetMask()
        void ClearSelected()
        void InvertMask() except +
        int NumAtomsInCommon(const _AtomMask&)
        void AddSelectedAtom(int i)
        void AddAtom(int)
        void AddAtoms(const vector [int]&)
        void AddAtomRange(int, int)
        void AddMaskAtPosition(const _AtomMask&, int)
        void PrintMaskAtoms(const char *)const 
        #int SetMaskString(const char *)
        int SetMaskString(const string&)
        void SetupIntMask(const char *, int, int)
        void SetupCharMask(const char *, int, int)
        bint AtomInCharMask(int)const 
        bint AtomsInCharMask(int, int)const 
        void SetNatom(int a)
        int ConvertToCharMask()
        int ConvertToIntMask()
        void MaskInfo()const 
        void BriefMaskInfo()const 
        #inline token_iterator begintoken()const 
        #inline token_iterator endtoken()const 

#ctypedef fused charstring:
#    char*
#    string

cdef class AtomMask:
    cdef _AtomMask* thisptr


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

cdef extern from "CpptrajState.h": 
    cdef cppclass _CpptrajState "CpptrajState":
        _CpptrajState()
        _TopologyList * PFL()
        _DataSetList * DSL()
        _DataFileList * DFL()
        void SetNoExitOnError()
        void SetNoProgress()
        void SetActionSilence(bint b)
        bint ExitOnError()const 
        bint EmptyState()const 
        int AddTrajin(_ArgList &, bint)
        int AddTrajin(const string&)
        int RunAnalyses()
        inline int AddTrajout "AddOutputTrajectory" (const _ArgList&)
        inline int AddTrajout "AddOutputTrajectory" (const string&)
        int AddReference(const string&, _ArgList &)
        inline int AddReference(const string&)
        inline int AddAction(DispatchAllocatorType, _ArgList &)
        inline int AddAnalysis(DispatchAllocatorType, _ArgList &)
        int WorldSize()
        int ListAll(_ArgList &)const 
        int SetListDebug(_ArgList &)
        int ClearList(_ArgList &)
        int RemoveDataSet(_ArgList &)
        int TrajLength(const string&, const vector[string]&)
        int Run()
        void MasterDataFileWrite()

cdef class CpptrajState:
    cdef _CpptrajState* thisptr
    cdef public TopologyList toplist
    cdef public DataFileList datafilelist
    cdef public DataSetList datasetlist

cdef extern from "ArgList.h": 
    cdef cppclass _ArgList "ArgList":
        _ArgList() 
        _ArgList(const char *)
        _ArgList(const string&)
        _ArgList(const string&, const char *)
        _ArgList(const _ArgList&)
        #_ArgList& operator =(const _ArgList&)
        const string& operator[](int) const 
        const vector[string]& List() const 
        vector[string].const_iterator begin() const 
        vector[string].const_iterator end() const 
        int Nargs() const 
        bint empty() const 
        const char * ArgLine() const 
        void ClearList() 
        int SetList(const string&, const char *)
        _ArgList RemainingArgs() 
        void AddArg(const string&)
        void MarkArg(int)
        bint CheckForMoreArgs() const 
        void PrintList() const 
        void PrintDebug() const 
        void RemoveFirstArg() 
        const char * Command() const 
        bint CommandIs(const char *) const 
        const string& GetStringNext() 
        const string& GetMaskNext() 
        const string& getNextTag() 
        bint ValidInteger(int) const 
        int IntegerAt(int) const 
        bint ValidDouble(int) const 
        int getNextInteger(int)
        double getNextDouble(double)
        const string& GetStringKey(const char *)
        int getKeyInt(const char *, int)
        double getKeyDouble(const char *, double)
        bint hasKey(const char *)
        bint Contains(const char *) const 

cdef class ArgList:
    cdef _ArgList* thisptr
