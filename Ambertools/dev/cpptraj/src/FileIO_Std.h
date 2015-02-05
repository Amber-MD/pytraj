#ifndef INC_FILEIO_STD_H
#define INC_FILEIO_STD_H
#include <cstdio> // For FILE
#include "FileIO.h"
// Class: FileIO_Std
/// File IO using CSTDLIB routines.
class FileIO_Std : public FileIO {
  public:
    FileIO_Std();
    ~FileIO_Std();
    int Open(const char *, const char *);    
    int Close();
    int Read(void *, size_t );
    int Write(const void *, size_t);
    int Flush()             { return fflush(fp_); }
    int Seek(off_t);
    int Rewind();  
    off_t Tell();  // NOTE: Tell may be unnecessary if only for size reporting.
    int Gets(char *, int );
    off_t Size(const char*) { return 0; }
    int SetSize(long int)   { return 0; }
  private:
    FILE *fp_;
    bool isStdout_;
};
#endif
