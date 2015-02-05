// CpptrajFile
#include <cstring>    // strlen 
#include <sys/stat.h> // stat
#include <cstdio>     // vsprintf
#include <cstdarg>    // va_X functions
#include "CpptrajFile.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // fileExists, tildeExpansion
// File Types
#include "FileIO_Std.h"
#ifdef HASGZ
#  include "FileIO_Gzip.h"
#endif
#ifdef MPI
#  include "FileIO_Mpi.h"
#endif
#ifdef HASBZ2
#  include "FileIO_Bzip2.h"
#endif

const char* CpptrajFile::FileTypeName[] = {
  "UNKNOWN_TYPE", "STANDARD", "GZIPFILE", "BZIP2FILE", "ZIPFILE", "MPIFILE"
};

// CONSTRUCTOR
CpptrajFile::CpptrajFile() :
  IO_(0),
  access_(READ),
  isDos_(0),
  uncompressed_size_(0UL),
  file_size_(0UL),
  compressType_(NO_COMPRESSION),
  debug_(0),
  isOpen_(false),
  fileType_(STANDARD)
{ }

// Copy Constructor
CpptrajFile::CpptrajFile(const CpptrajFile &rhs) :
  IO_(0),
  access_(rhs.access_),
  isDos_(rhs.isDos_),
  uncompressed_size_(rhs.uncompressed_size_),
  file_size_(rhs.file_size_),
  compressType_(rhs.compressType_),
  debug_(rhs.debug_),
  // Even if file is open, copy it closed?
  isOpen_(false),
  fileType_(rhs.fileType_),
  fname_(rhs.fname_)
{
  // Set up the IO object
  // NOTE: Should probably throw an exception if this fails.
  if (rhs.IO_ != 0) 
    IO_ = SetupFileIO( fileType_ );
}

// Assignment
CpptrajFile &CpptrajFile::operator=(const CpptrajFile &rhs) {
  // Self-assignment
  if (this == &rhs) return *this;
  // Deallocate
  CloseFile();
  if (IO_ != 0) delete IO_;
  // Allocate and copy
  debug_ = rhs.debug_;
  access_ = rhs.access_;
  uncompressed_size_ = rhs.uncompressed_size_;
  file_size_ = rhs.file_size_;
  fileType_ = rhs.fileType_;
  fname_ = rhs.fname_;
  compressType_ = rhs.compressType_;
  isDos_ = rhs.isDos_;
  // Set up the IO object
  // NOTE: Should probably throw an exception if this fails.
  if (rhs.IO_ != 0) 
    IO_ = SetupFileIO( fileType_ );
  else
    IO_ = 0;
  return *this;
}

// DESTRUCTOR
CpptrajFile::~CpptrajFile() {
   //fprintf(stderr,"CPPTRAJFILE DESTRUCTOR\n");
   CloseFile();
   if (IO_ != 0) delete IO_;
}

// CpptrajFile::IsCompressed()
bool CpptrajFile::IsCompressed() {
  if (compressType_ == NO_COMPRESSION) return false;
  return true;
}

int CpptrajFile::OpenFile() {
  return OpenFile( access_ );
}

// CpptrajFile::OpenFile()
/** Open the file. If already open, reopen.  */
int CpptrajFile::OpenFile(AccessType accessIn) {
  int err;
  if (isOpen_) CloseFile();

  switch (accessIn) {
    case READ:
      if (fname_.empty()) {
        mprinterr("Error: CpptrajFile: Filename is null.\n");
        return 1;
      }
      if ( IO_->Open(Filename().full(), "rb")  ) { 
        rprintf("Could not open %s for reading.\n", Filename().full());
        return 1;
      }
      if (debug_>0) rprintf("Opened %s for reading.\n", Filename().full());
      break;
    case APPEND:
      if (fname_.empty()) {
        mprinterr("Error: CpptrajFile: Filename is null.\n");
        return 1;
      }
      if ( IO_->Open(Filename().full(), "ab") ) {
        rprintf("Could not open %s for appending.\n", Filename().full());
        return 1;
      }
      if (debug_>0) rprintf("Opened %s for appending.\n", Filename().full());
      break;
    case WRITE:
      err = 0;
      if ( fname_.empty() )
        err = IO_->Open(0, "wb");
      else
        err = IO_->Open(Filename().full(), "wb");
      if ( err != 0 ) { 
        rprintf("Could not open %s for writing.\n", Filename().full());
        return 1;
      }
      if (debug_>0) rprintf("Opened %s for writing.\n", Filename().full());
      break;
    case UPDATE:
      if (fname_.empty()) {
        mprinterr("Error: CpptrajFile: Filename is null.\n");
        return 1;
      }
      if ( IO_->Open(Filename().full(), "r+b") ) {
        rprintf("Could not open %s for updating.\n", Filename().full());
        return 1;
      }
      if (debug_>0) rprintf("Opened %s for updating.\n", Filename().full());
      break;
  }
      
  isOpen_ = true;
  return 0;
}

// CpptrajFile::CloseFile()
void CpptrajFile::CloseFile() {
  if (isOpen_) {
    IO_->Close();
    if (debug_>0) rprintf("Closed %s.\n", Filename().full());
    isOpen_=false;
  }
}

// -----------------------------------------------------------------------------
// CpptrajFile::Printf()
/** Take the formatted string and write it to file using Write.
  */
void CpptrajFile::Printf(const char *format, ...) {
  va_list args;
  va_start(args, format);
  vsprintf(linebuffer_,format,args);
  IO_->Write(linebuffer_, strlen(linebuffer_));
  va_end(args);
}

std::string CpptrajFile::GetLine() {
  if (IO_->Gets(linebuffer_, BUF_SIZE) != 0) {
    //mprinterr("Error: Getting line from %s\n", Filename().full());
    return std::string();
  }
  return std::string(linebuffer_);
}

const char* CpptrajFile::NextLine() {
  if (IO_->Gets(linebuffer_, BUF_SIZE) != 0) {
    //mprinterr("Error: Reading line from %s\n", Filename().full());
    return 0;
  }
  return linebuffer_;
}

// -----------------------------------------------------------------------------
// CpptrajFile::UncompressedSize()
off_t CpptrajFile::UncompressedSize() {
  if (compressType_ == NO_COMPRESSION)
    return file_size_;
  else
    return uncompressed_size_;
}

// CpptrajFile::Reset()
/** Close file if open, reset all file information.
  */
void CpptrajFile::Reset() {
  CloseFile();
  if (IO_!=0) delete IO_;
  IO_ = 0;
  fname_.clear();
  isOpen_ = false;
  uncompressed_size_ = 0UL;
  compressType_ = NO_COMPRESSION;
  isDos_ = 0;
}

// CpptrajFile::OpenRead()
int CpptrajFile::OpenRead(std::string const& nameIn) {
  if (SetupRead( nameIn, debug_ )) return 1; 
  return OpenFile();
}

// CpptrajFile::SetupRead()
/** Set up file for reading. Will autodetect the type.
  * \return 0 on success, 1 on error.
  */
int CpptrajFile::SetupRead(std::string const& nameIn, int debugIn) {
  // null filename not allowed
  if (nameIn.empty()) {
    mprinterr("Internal Error: No filename specified for READ.\n");
    return 1;
  }
  // Check if file exists. If not, fail silently
  if (!fileExists( nameIn )) return 1;
  // Clear file, set debug level
  Reset();
  debug_ = debugIn;
  access_ = READ;
  fileType_ = UNKNOWN_TYPE;
  if (debug_>0)
    mprintf("CpptrajFile: Setting up %s for READ.\n", nameIn.c_str());
  // Perform tilde-expansion
  std::string expandedName = tildeExpansion( nameIn );
  if (expandedName.empty()) {
    mprinterr("Interal Error: CpptrajFile: Tilde-expansion failed.\n");
    return 1;
  }
  // Determine file type. This sets up IO and determines compression. 
  if (ID_Type( expandedName.c_str() )) return 1;
  // Set up filename; sets base filename and extensions
  fname_.SetFileName( expandedName, IsCompressed() );
  if (debug_>0)
    rprintf("\t[%s] is type %s with access READ\n", Filename().full(), FileTypeName[fileType_]);
  return 0;
}

// CpptrajFile::OpenWriteNumbered()
// NOTE: File MUST be previously set up. Primarily for use with traj files.
int CpptrajFile::OpenWriteNumbered(int numIn) {
  std::string newName = NumberFilename( Filename().Full(), numIn );
  if (IO_->Open( newName.c_str(), "wb")) return 1;
  isOpen_ = true;
  return 0;
}

// CpptrajFile::OpenWrite()
int CpptrajFile::OpenWrite(std::string const& nameIn) {
  if (SetupWrite(nameIn, debug_)) return 1;
  return OpenFile();
}

// CpptrajFile::OpenEnsembleWrite()
int CpptrajFile::OpenEnsembleWrite(std::string const& nameIn, int ensembleNum) {
  if (!nameIn.empty() && ensembleNum > -1) {
    if (SetupWrite( NumberFilename(nameIn, ensembleNum), debug_)) return 1;
  } else {
    if (SetupWrite( nameIn,                              debug_)) return 1;
  }
  return OpenFile();
}

// CpptrajFile::SetupWrite()
int CpptrajFile::SetupWrite(std::string const& nameIn, int debugIn) {
  return SetupWrite(nameIn, UNKNOWN_TYPE, debugIn);
}

// CpptrajFile::SetupWrite()
/** Set up file for writing with the given type. If no filename is given 
  * this indicates STDOUT. If no type is specified attempt to detect
  * from the compression extension.
  * \return 0 on success, 1 on error.
  */
int CpptrajFile::SetupWrite(std::string const& filenameIn, FileType typeIn, int debugIn) 
{
  // Clear file, set debug level
  Reset();
  debug_ = debugIn;
  access_ = WRITE;
  if (debug_>0)
    mprintf("CpptrajFile: Setting up %s for WRITE.\n",filenameIn.c_str());
  // Set up filename; sets base filename and extension
  fname_.SetFileName(filenameIn);
  // If file type is not specified, try to determine from filename extension
  if (typeIn == UNKNOWN_TYPE) {
    if (fname_.Compress() == ".gz")
      fileType_ = GZIPFILE;
    else if (fname_.Compress() == ".bz2")
      fileType_ = BZIP2FILE;
    else
      fileType_ = STANDARD;
  }
  // Setup IO based on type.
  IO_ = SetupFileIO( fileType_ );
  if (IO_ == 0) return 1;
  if (debug_>0)
    rprintf("\t[%s] is type %s with access WRITE\n", Filename().full(), FileTypeName[fileType_]);
  return 0;
}

// CpptrajFile::OpenAppend()
int CpptrajFile::OpenAppend(std::string const& nameIn) {
  if (nameIn.empty()) {
    if (SetupWrite(nameIn, debug_)) return 1;
  } else {
    if (SetupAppend(nameIn, debug_)) return 1;
  }
  return OpenFile();
}

// CpptrajFile::OpenEnsembleAppend()
int CpptrajFile::OpenEnsembleAppend(std::string const& nameIn, int ensembleNum) {
  if (nameIn.empty())
    return OpenEnsembleWrite( nameIn, ensembleNum );
  else {
    if (ensembleNum > -1) {
      if (SetupAppend( NumberFilename(nameIn, ensembleNum), debug_)) return 1;
    } else {
      if (SetupAppend( nameIn,                              debug_)) return 1;
    }
  }
  return OpenFile();
}

// CpptrajFile::SetupAppend()
/** Set up the file for appending. Will first set up for read to determine
  * the type and format.
  * \return 0 on success, 1 on error.
  */
int CpptrajFile::SetupAppend(std::string const& filenameIn, int debugIn) {
  // Make append to null an error
  if (filenameIn.empty()) {
    mprinterr("Error: SetupAppend(): No filename specified\n");
    return 1;
  }
  // NOTE: File will be cleared and debug set by either SetupRead/SetupWrite
  if (fileExists(filenameIn)) {
    // If file exists, first set up for read to determine type and format.
    if (SetupRead(filenameIn, debugIn)!=0) return 1;
    access_ = APPEND;
  } else {
    // File does not exist, just set up for write.
    if (SetupWrite(filenameIn, debugIn)!=0) return 1;
    if (debug_>0)
      mprintf("Warning: %s does not exist, changed access from APPEND to WRITE.\n",
              Filename().full());
  }
  // Appending and compression not supported.
  if (IsCompressed()) {
    mprinterr("Error: Appending to compressed files is not supported.\n");
    return 1;
  }
  return 0;
}

// -----------------------------------------------------------------------------
// CpptrajFile::SetupFileIO()
/** Set up the IO based on given file type. */
FileIO* CpptrajFile::SetupFileIO(FileType typeIn) {
  switch (typeIn) {
    case STANDARD  : return (new FileIO_Std());
    case GZIPFILE  : 
#ifdef HASGZ
      return new FileIO_Gzip(); 
#else
      mprinterr("Error: Compiled without Gzip support. Recompile with -DHASGZ\n");
      return 0;
#endif
      break;
    case BZIP2FILE :
#ifdef HASBZ2 
      return (new FileIO_Bzip2());
#else
      mprinterr("Error: Compiled without Bzip2 support. Recompile with -DHASBZ2\n");
      return 0;
#endif
    break;
    case MPIFILE   : 
#ifdef MPI
      return (new FileIO_Mpi());
#else
      mprinterr("Error: Compiled without MPI support. Recompile with -DMPI\n");
      return 0;
#endif
      break;
    //case ZIPFILE   : return (new ZipFile()); break;
    default : 
      mprinterr("Error: Unrecognized file type.\n");
  }
  return 0;
}

// CpptrajFile::ID_Type() 
/** Attempt to identify the file type for filenameIn. Also set file_size,
  * uncompressed_size, and compressType.
  */
int CpptrajFile::ID_Type(const char* filenameIn) {
  if (filenameIn == 0) return 1;
  // Get basic file information
  struct stat frame_stat;
  if (stat(filenameIn, &frame_stat) == -1) {
    mprinterr( "Error: Could not find file status for %s\n", filenameIn);
    if (debug_>0) 
      perror("     Error from stat: ");
    return 1;
  }
  file_size_ = frame_stat.st_size;
  // Start off every file as a standard file
  fileType_ = STANDARD;
  IO_ = new FileIO_Std();
  // ID by magic number - open for binary read access
  if ( IO_->Open(filenameIn, "rb") ) { 
    mprintf("Could not open %s for hex signature read.\n", filenameIn);
    return 1;
  }
  // Read first 3 bytes
  unsigned char magic[3];
  magic[0] = 0; 
  magic[1] = 0; 
  magic[2] = 0;
  IO_->Read(magic, 3);
  IO_->Close();
  if (debug_>0) mprintf("\t    Hex sig: %x %x %x", magic[0],magic[1],magic[2]);
  // Check compression
  if ((magic[0]==0x1f) && (magic[1]==0x8b) && (magic[2]==0x8)) {
    if (debug_>0) mprintf(", Gzip file.\n");
    compressType_ = GZIP;
    fileType_ = GZIPFILE;
  } else if ((magic[0]==0x42) && (magic[1]==0x5a) && (magic[2]==0x68)) {
    if (debug_>0) mprintf(", Bzip2 file.\n");
    compressType_ = BZIP2;
    fileType_ = BZIP2FILE;
  } else if ((magic[0]==0x50) && (magic[1]==0x4b) && (magic[2]==0x3)) {
    if (debug_>0) mprintf(", Zip file.\n");
    compressType_ = ZIP;
    fileType_ = ZIPFILE;
  } else {
    if (debug_>0) mprintf(", No compression.\n");
  }
  // Assign the appropriate IO type based on the file type
  delete (FileIO_Std*) IO_;
  IO_ = SetupFileIO( fileType_ );
  if (IO_ == 0) return 1;

  // If the file is compressed, get the uncompressed size
  // For standard files this just returns 0UL
  // Standard file size is in the frame_stat struct
  uncompressed_size_ = IO_->Size(filenameIn);

  // Additional file characteristics
  linebuffer_[0]='\0';
  if (IO_->Open(filenameIn, "rb")!=0) return 1; 
  IO_->Gets(linebuffer_,BUF_SIZE);
  IO_->Close();

  // Check for terminal CR before newline, indicates DOS file
  size_t i = strlen(linebuffer_);
  if ( i>1 ) {
    if (linebuffer_[ i - 2 ] == '\r') {
      if (debug_>0) mprintf("  [DOS]");
      isDos_ = 1;
    }
  }
  return 0;
}
