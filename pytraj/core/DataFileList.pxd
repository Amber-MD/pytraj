# distutils: language = c++
from libcpp.vector cimport vector
from libcpp.string cimport string
from pytraj.DataFile cimport _DataFile, DataFile
from pytraj.datasets.DataSet cimport _DataSet, DataSet
from pytraj.ArgList cimport _ArgList, ArgList


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
