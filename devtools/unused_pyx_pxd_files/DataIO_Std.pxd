# distutils: language = c++
from DataIO cimport *


cdef extern from "DataIO_Std.h": 
    cdef cppclass _DataIO_Std "DataIO_Std":
        _DataIO_Std() 
        _BaseIOtype * Alloc() 
        void ReadHelp() 
        void WriteHelp() 
        int ReadData(const string&, _ArgList &, _DataSetList &, const string&)
        int processWriteArgs(_ArgList &)
        int WriteData(const string&, const _DataSetList&)
        int WriteData2D(const string&, const _DataSetList&)
        int WriteData3D(const string&, const _DataSetList&)
        bint ID_DataFormat(_CpptrajFile &)
