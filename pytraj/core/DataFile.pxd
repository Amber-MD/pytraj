# distutils: language = c++
from __future__ import absolute_import
from libcpp.string cimport string
from ..ArgList cimport _ArgList, ArgList
from .cpptraj_core cimport _FileName, FileName

# for some reasons, I need to use absolute path here
from ..datasets.base cimport _DataSet, DataSet
from ..datasets.DataSetList cimport _DataSetList, DataSetList


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
