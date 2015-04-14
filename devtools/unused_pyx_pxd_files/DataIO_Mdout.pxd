# distutils: language = c++
from DataIO cimport *


cdef extern from "DataIO_Mdout.h": 
    cdef cppclass _DataIO_Mdout "DataIO_Mdout":
        _DataIO_Mdout() 
        _BaseIOtype * Alloc() 
        void ReadHelp() 
        int ReadData(const string&, _ArgList&, _DataSetList&, const string&)
        int processWriteArgs(_ArgList&)
        int WriteData(const string&, const _DataSetList&)
        int WriteData2D(const string&, const _DataSetList&)
        int WriteData3D(const string&, const _DataSetList&)
        bint ID_DataFormat(_CpptrajFile&)

cdef extern from "DataIO_Mdout.cpp":
    inline int EOF_ERROR()
    #inline int EOF_ERROR()
