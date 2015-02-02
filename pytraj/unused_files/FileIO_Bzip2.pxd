# distutils: language = c++
from posix.unistd cimport off_t
#from bzlib cimport *
from FileIO cimport *


cdef extern from "FileIO_Bzip2.h": 
    cdef cppclass _FileIO_Bzip2 "FileIO_Bzip2":
        _FileIO_Bzip2() 
        #~_FileIO_Bzip2() 
        int Open(const char *, const char *)
        int Close() 
        off_t Size(const char *)
        int Read(void *, size_t)
        int Write(const void *, size_t)
        int Flush() 
        int Seek(off_t)
        int Rewind() 
        off_t Tell() 
        int Gets(char *, int)
        int SetSize(long int)
