# distutils: language = c++
from libcpp.string cimport string


cdef extern from "ReadLine.h": 
    cdef cppclass _ReadLine "ReadLine":
        _ReadLine() 
        int GetInput() 
        bint YesNoPrompt(const char *)
        const char * c_str() const 
        const string& operator *() const 
        bint empty() const 
