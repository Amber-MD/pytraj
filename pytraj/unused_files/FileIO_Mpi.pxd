# distutils: language = c++
from posix.unistd cimport off_t
from FileIO cimport *
from MpiRoutines cimport *


cdef extern from "FileIO_Mpi.h": 
    cdef cppclass _FileIO_Mpi "FileIO_Mpi":
        _FileIO_Mpi() 
        #~_FileIO_Mpi() 
        int Open(const char *, const char *)
        int Close() 
        int Read(void *, size_t)
        int Write(const void *, size_t)
        int Flush() 
        int Seek(off_t)
        int Rewind() 
        off_t Tell() 
        int Gets(char *, int)
        int SetSize(long int)
        off_t Size(const char *)
