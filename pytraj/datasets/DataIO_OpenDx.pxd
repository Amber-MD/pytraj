# distutils: language = c++
from DataIO cimport *


cdef extern from "DataIO_OpenDx.h": 
    cdef cppclass _DataIO_OpenDx "DataIO_OpenDx":
        _DataIO_OpenDx() : DataIO(false, false, true)
        _BaseIOtype * Alloc() 
        int ReadData(const string&, _ArgList&, _DataSetList&, const string&)
        int processWriteArgs(_ArgList&)
        int WriteData(const string&, const _DataSetList&)
        int WriteData2D(const string&, const _DataSetList&)
        int WriteData3D(const string&, const _DataSetList&)
        bint ID_DataFormat(_CpptrajFile&)
