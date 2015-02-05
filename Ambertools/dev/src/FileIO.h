#ifndef INC_FILEIO_H
#define INC_FILEIO_H
#include <sys/types.h> // For off_t
// Class: FileIO 
/** Abstract base class for handling basic file IO. File types will inherit 
  * this class and replace basic IO functions (open, close, read etc) with
  * their own.
  * NOTES:
  *   File pointers for each file type will be declared in their own class.
  *   Use off_t for seeking for better compatibility with large files.
  *   Thus far it has only been necessary to seek using SEEK_SET, so that
  *   argument will be removed from Seek() for now.
  */ 
class FileIO {
  public:
    virtual ~FileIO() {}
    /// Open the file with given name and mode.
    virtual int Open(const char *, const char *) = 0;
    /// Close the file.
    virtual int Close()                          = 0;
    /// Read bytes from a file.
    /** \return number of bytes read, -1 on error.
      */
    virtual int Read(void *, size_t)             = 0;
    /// Write bytes to a file.
    /** \return 0 on success, 1 on failure.
      */
    virtual int Write(const void *, size_t)      = 0;
    /// Flush file to disk.
    virtual int Flush()                          = 0;
    /// Seek to specified position in file.
    virtual int Seek(off_t)                      = 0;
    /// Reset file pointer to beginning.
    virtual int Rewind()                         = 0;
    /// Return current file position.
    virtual off_t Tell()                         = 0;
    /// Get a line from the file.
    virtual int Gets(char *, int)                = 0;
    /// Return uncompressed file size; only needed for compressed files.
    virtual off_t Size(const char *)             = 0;
    /// Set expected file size, only needed for writes with MPI files.
    virtual int SetSize(long int)                = 0;
};
#endif
