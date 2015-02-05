#ifndef INC_FILEIO_MPI_H
#define INC_FILEIO_MPI_H
#include "FileIO.h" 
#include "MpiRoutines.h"
// Class: FileIO_Mpi
/// MPI file IO, wrappers for the MPI routines in MpiRoutines.h
class FileIO_Mpi : public FileIO {
  public:
    FileIO_Mpi(); 
    ~FileIO_Mpi(); 

    int Open(const char *, const char *);    
    int Close();
    int Read(void *, size_t );
    int Write(const void *, size_t);
    int Flush();
    int Seek(off_t);
    int Rewind();  
    off_t Tell();
    int Gets(char *, int );
    int SetSize(long int);
    off_t Size(const char*) { return 0; }
  private:
    parallelType pfile_; 
};
#endif
