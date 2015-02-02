# distutils: language = c++
from DataIO cimport *


cdef extern from "DataIO_Xplor.h": 
    cdef cppclass _DataIO_Xplor "DataIO_Xplor":
        _DataIO_Xplor() : DataIO(false, false, true)
        _BaseIOtype * Alloc() 
        int ReadData(const string&, _ArgList &, _DataSetList &, const string&)
        int processWriteArgs(_ArgList &)
        int WriteData(const string&, const _DataSetList&)
        int WriteData2D(const string&, const _DataSetList&)
        int WriteData3D(const string&, const _DataSetList&)
        bint ID_DataFormat(_CpptrajFile &)
