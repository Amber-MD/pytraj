#ifndef INC_CPPTRAJFILE_H
#define INC_CPPTRAJFILE_H
#include "FileName.h" 
#include "FileIO.h"
// Class: CpptrajFile
/// Class to abstract handling of basic file routines.
class CpptrajFile {
  public:
    enum AccessType   { READ, WRITE, APPEND, UPDATE };
    enum CompressType { NO_COMPRESSION, GZIP, BZIP2, ZIP };
    enum FileType { UNKNOWN_TYPE, STANDARD, GZIPFILE, BZIP2FILE, ZIPFILE, MPIFILE };

    CpptrajFile();
    virtual ~CpptrajFile(); // Virtual since class is inherited
    CpptrajFile(const CpptrajFile&);
    CpptrajFile &operator=(const CpptrajFile &);
    /// Set up and open file for reading.
    int OpenRead(std::string const&);
    /// Prepare file for reading. 
    int SetupRead(std::string const&, int);
    /// Open file previously set up for write with given numeric suffix.
    int OpenWriteNumbered(int);
    /// Set up and open file for writing.
    int OpenWrite(std::string const&);
    /// Open file for writing in Actions; if ensemble, append numeric suffix to filename.
    int OpenEnsembleWrite(std::string const&, int);
    /// Prepare file for writing.
    int SetupWrite(std::string const&, int);
    /// Prepare file of given type for writing
    int SetupWrite(std::string const&, FileType, int);
    /// Set up and open file for appending.
    int OpenAppend(std::string const&);
    /// Open file for append in Actions; if ensemble, append numeric suffix to filename.
    int OpenEnsembleAppend(std::string const&, int);
    /// Prepare file for appending. 
    int SetupAppend(std::string const&, int);
    /// Open file.
    int OpenFile();
    /// Open file with given access.
    int OpenFile(AccessType);
    /// Close file.
    void CloseFile();
    /// Printf using the files Write routine.
    void Printf(const char*, ...);
    /// Get next line as a string
    std::string GetLine();
    /// Get next line and return pointer to raw buffer
    const char* NextLine();
    /// \return the access file is currently set up for.
    AccessType Access()         const { return access_;               }
    /// \return the compression type
    CompressType Compression()  const { return compressType_;         }
    /// \return true if the file is open
    bool IsOpen()               const { return isOpen_;               }
    /// \return file name class.
    const FileName& Filename()  const { return fname_;                }
    /// \return 1 if the file contains carriage returns in addition to newlines
    int IsDos()                 const { return isDos_;                }
    /// \return File size
    off_t FileSize()            const { return file_size_;            }
    /// \return true if the file is compressed.
    bool IsCompressed();
    /// \return uncompressed file size (just size if file is not compressed).
    off_t UncompressedSize();
    int Gets(char* buf, int num)           { return IO_->Gets(buf, num);  }
    int Write(const void* buf, size_t num) { return IO_->Write(buf, num); }
    int Read(void* buf, size_t num)        { return IO_->Read(buf, num);  }
    int Seek(off_t offset)                 { return IO_->Seek(offset);    }
    int Rewind()                           { return IO_->Rewind();        }
    int Flush()                            { return IO_->Flush();         }
    off_t Tell()                           { return IO_->Tell();          }
  protected: // Protected for PDBfile/Mol2File
    static const size_t BUF_SIZE = 1024;
    char linebuffer_[BUF_SIZE]; ///< Used in Printf/GetLine functions
  private:
    static const char* FileTypeName[];
    FileIO* IO_;                ///< The interface to basic IO operations.
    AccessType access_;         ///< Access (Read, write, append)
    int isDos_;                 ///< 1 if CR present, need to count them as newlines
    off_t uncompressed_size_;   ///< If compressed, uncompressed file size
    off_t file_size_;           ///< Actual file size
    CompressType compressType_; ///< Type of compression
    int debug_;                 ///< Debug level
    bool isOpen_;               ///< If true, file is open and ready for IO.
    FileType fileType_;         ///< File type (determines IO)
    FileName fname_;            ///< Holds full and base file name + any extensions.

    void Reset();
    static FileIO* SetupFileIO(FileType);
    int ID_Type(const char*);
};
#endif
