# distutils: language = c++
from libcpp.string cimport string


cdef extern from "DataSet.h": 
    ctypedef enum DataType "DataSet::DataType":
        pass
    ctypedef enum scalarMode "DataSet::scalarMode":
        pass
    ctypedef enum scalarType "DataSet::scalarType":
        pass
    cdef cppclass _DataSet "DataSet":
        _DataSet() 
        _DataSet(DataType, int, int, int)
        _DataSet(const _DataSet&)
        #_DataSet& operator =(const _DataSet&)
        void SetWidth(int)
        void SetPrecision(int, int)
        int SetupSet(const string&, int, const string&, int)
        void SetLegend(const string& lIn)
        void SetScalar(scalarMode mIn)
        inline void SetScalar(scalarMode, scalarType)
        int SetDataSetFormat(bint)
        bint Matches(const string&, int, const string&) const 
        void ScalarDescription() const 
        bint Empty() const 
        const string& Legend() const 
        const string& Name() const 
        int Idx() const 
        const string& Aspect() const 
        int ColumnWidth() const 
        DataType Type() const 
        scalarMode ScalarMode() const 
        scalarType ScalarType() const 
        size_t Ndim() const 
        inline bint operator< (const _DataSet&) const 
        const char * DataFormat() const 
        scalarMode ModeFromKeyword(const string&)
        scalarType TypeFromKeyword(const string&, scalarMode&)
        scalarType TypeFromKeyword(const string&, const scalarMode&)
        size_t Size()

        # virtual
        #void Add( size_t, const void*  )

cdef class DataSet:
    cdef _DataSet* baseptr0
    cdef public object _base
