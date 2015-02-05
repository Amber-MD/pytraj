#ifndef INC_FILEIO_BZIP2_H
#define INC_FILEIO_BZIP2_H
#ifdef HASBZ2
// NOTE: bzlib.h has stdio.h. Does it matter that its not cstdio? 
#include "bzlib.h"
#include "FileIO.h"
// Class: FileIO_Bzip2
/// Bzip2 file IO.
class FileIO_Bzip2 : public FileIO {
  public:
    FileIO_Bzip2(); 
    ~FileIO_Bzip2(); 
    int Open(const char *, const char *);    
    int Close();
    off_t Size(const char *);
    int Read(void *, size_t );
    int Write(const void *, size_t);
    int Flush() { return 0; } // bzflush is not part of the bzip standard
    int Seek(off_t);
    int Rewind();  
    off_t Tell();  // NOTE: Tell may be unnecessary if only for size reporting.
    int Gets(char *, int );
    int SetSize(long int) { return 0; }
  private:
    bool isBzread_;
    FILE *fp_;
    BZFILE *infile_;
    int err_;
    char *bzfilename_;
    char *bzmode_;
    off_t position_;

    const char *BZerror();
};
#endif
#endif
