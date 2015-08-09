# distutils: language = c++
from DataIO cimport *


cdef extern from "DataIO_Grace.h": 
    cdef cppclass _DataIO_Grace "DataIO_Grace":
        _DataIO_Grace() : DataIO(true, false, false ), isInverted_(false)
        _BaseIOtype * Alloc() 
        void WriteHelp() 
        int ReadData(const string&, _ArgList &, _DataSetList &, const string&)
        int processWriteArgs(_ArgList &)
        int WriteData(const string&, const _DataSetList&)
        int WriteData2D(const string&, const _DataSetList&)
        int WriteData3D(const string&, const _DataSetList&)
        bint ID_DataFormat(_CpptrajFile &)
