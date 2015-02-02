# distutil: language = c++

from libcpp.string cimport string

cdef extern from "NameType.h":
    cdef cppclass _NameType "NameType":
        _NameType() 
        _NameType(const _NameType&)
        _NameType(const char *)
        _NameType(const string&)
        #_NameType& operator =(const _NameType&)
        void ToBuffer(char *) const 
        bint Match(const _NameType&) const 
        bint operator ==(const _NameType&) const 
        bint operator ==(const char *) const 
        #bint opr_ne "operator !="(const _NameType&) const 
        bint operator !=(const _NameType&) const 
        bint operator !=(const char *) const 
        const char* opr_star "operator*" () const 
        char opr_idx "operator[]"(int) const 
        string Truncated() const 
        void ReplaceAsterisk() 


cdef class NameType:
        cdef _NameType* thisptr
