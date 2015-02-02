# distutils: language = c++

from libcpp.string cimport string

cdef extern from "FileName.h":
    cdef cppclass _FileName "FileName":
        _FileName()
        _FileName(_FileName)
        int SetFileName(string)
        int SetFileNameWithExpansion(string)
        int SetFileName(string, bool)
        void clear()
        bint MatchFullOrBase(string)
        string Full()
        string Base()
        char * full()
        char * base()
        string Ext()
        string Compress()
        string DirPrefix()
        bint empty()

cdef class FileName:
    cdef _FileName* thisptr
