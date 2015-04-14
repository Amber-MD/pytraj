# distutil: language = c++

from libcpp.string cimport string
cdef extern from "BufferedLine.h":
    cdef cppclass _BufferedLine "BufferedLine":
        _BufferedLine()
        const char* Line()
        int TokenizeLine(const char*)
        const char* NextToken()
        inline const char* Token(int)
        int OpenFileRead(const string& fname)
        int LineNumber()
        const char* Buffer()
        
cdef class BufferedLine:
    cdef _BufferedLine* thisptr
