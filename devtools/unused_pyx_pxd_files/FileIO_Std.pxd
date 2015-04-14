# distutil: language = c++
from posix.unistd cimport off_t
from FileIO cimport *
cdef extern from "FileIO_Std.h":
    cdef cppclass _FileIO_Std "FileIO_Std":
        _FileIO_Std()
        int Open(const char*, const char*)
        int Close()
        int Read(void*, size_t)
        int Write(const void*, size_t)
        int Flush()
        int Seek(off_t)
        int Rewind()
        off_t Tell()
        int Gets(char*, int)
        off_t Size(const char*)
        int SetSize(long int)
