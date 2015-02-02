# distutils: language = c++
from libcpp.string cimport string


cdef extern from "Range.h": 
    cdef cppclass _Range "Range":
        _Range() 
        _Range(const string&)
        _Range(const string&, int)
        _Range(const _Range&)
        #_Range& operator =(const _Range&)
        #const_iterator begin() const 
        #const_iterator end() const 
        bint Empty() const 
        int Size() const 
        int Set_Range(const string&)
        int Set_Range(int, int)
        const char * RangeArg() const 
        void Print_Range(const char *, int) const 
        void ShiftBy(int)
        void AddTo_Range(int)
        void RemoveFrom_Range(int)
        bint In_Range(int) const 
