# distutils: language = c++
from libcpp.string cimport string
from ..datasets.DataSet cimport _DataSet, DataSet
from ..ArgList cimport _ArgList, ArgList
from ..datasets.DataSetList cimport _DataSetList, DataSetList
from .FileName cimport _FileName, FileName


cdef extern from "DataFile.h": 
    # DataFile.h
    ctypedef enum DataFormatType "DataFile::DataFormatType":
        pass
    cdef cppclass _DataFile "DataFile":
        _DataFile() 
        #~_DataFile() 
        @staticmethod
        void WriteHelp() 
        @staticmethod
        void ReadOptions() 
        @staticmethod
        void WriteOptions() 
        @staticmethod
        DataFormatType GetFormatFromArg(_ArgList& a)
        #@staticmethod
        const char * FormatString(DataFormatType t)
        const char * FormatString() const 
        void SetDebug(int)
        void Set_DataFilePrecision(int, int)
        int ReadDataIn(const string&, const _ArgList&, _DataSetList&)
        int SetupDatafile(const string&, _ArgList&, int)
        void SetDataFilePrecision(int, int)
        int AddSet(_DataSet *)
        int RemoveSet(_DataSet *)
        int ProcessArgs(_ArgList&)
        int ProcessArgs(const string&)
        void WriteData() 
        void DataSetNames() const 
        const _FileName& _DataFilename() const 
        void SetDFLwrite(bint fIn)
        bint DFLwrite() const 
        DataFormatType Type() const 


cdef class DataFile:
    cdef _DataFile* thisptr
    cdef bint py_free_mem
