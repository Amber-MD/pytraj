# distutils: language = c++
from DataIO cimport *


cdef extern from "DataIO_Gnuplot.h": 
    cdef cppclass _DataIO_Gnuplot "DataIO_Gnuplot":
        _DataIO_Gnuplot() 
        _BaseIOtype * Alloc() 
        void WriteHelp() 
        int ReadData(const string&, _ArgList &, _DataSetList &, const string&)
        int processWriteArgs(_ArgList &)
        int WriteData(const string&, const _DataSetList&)
        int WriteData2D(const string&, const _DataSetList&)
        int WriteData3D(const string&, const _DataSetList&)
        bint ID_DataFormat(_CpptrajFile &)
