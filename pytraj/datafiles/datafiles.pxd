# distutils: language = c++
from __future__ import absolute_import
from libcpp.string cimport string
from libcpp.vector cimport vector

from ..datasets.base cimport _DataSet, DataSet
from ..core.cpptraj_core cimport _ArgList, ArgList, _ArgList, ArgList, _FileName, FileName

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


cdef extern from "DataFileList.h": 
    cdef cppclass _DataFileList "DataFileList":
        _DataFileList() 
        #~_DataFileList() 
        void Clear() 
        _DataFile * RemoveDataFile(_DataFile *)
        void RemoveDataSet(_DataSet *)
        void SetDebug(int)
        # this method is for MPI
        void SetEnsembleMode(int mIn)
        _DataFile * GetDataFile(const string&) const 
        _DataFile * AddDataFile(const string&, _ArgList&)
        _DataFile * AddDataFile(const string&)
        _DataFile * AddSetToFile(const string&, _DataSet *)
        void List() const 
        void WriteAllDF() 
        void ResetWriteStatus() 
        int ProcessDataFileArgs(_ArgList&)


cdef class DataFileList:
    cdef _DataFileList* thisptr
    cdef bint py_free_mem
