# distutil: language = c++
from posix.unistd cimport off_t

cdef extern from "FileIO_Gzip.h":
    cdef cppclass _FileIO_Gzip "FileIO_Gzip":
        _FileIO_Gzip()
        int Open(const char*, const char*)
        int Close()
        off_t Size(const char*)
        pass 
