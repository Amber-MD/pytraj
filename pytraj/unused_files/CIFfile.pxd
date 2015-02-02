# distutils: language = c++
from Atom cimport *
from BufferedLine cimport *


cdef extern from "CIFfile.h": 
    cdef cppclass _CIFfile "CIFfile":
        DataBlock() 
        const string& Header() const 
        bint empty() const 
        int AddHeader(const string&)
        int AddSerialDataRecord(const char *, _BufferedLine &)
        int AddLoopColumn(const char *)
        int AddLoopData(const char *, _BufferedLine &)
        void ListData() const 
        int ColumnIndex(const string&)const 
        string Data(const string&)const 
        #data_it begin() const 
        #data_it end() const 
