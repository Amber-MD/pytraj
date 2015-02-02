// FileIO_Bzip2: Bzip2 file operations
#ifdef HASBZ2
#include <cstring>
#include <cstdlib>
#include "FileIO_Bzip2.h" // FileIO.h, cstdio, bzlib.h
#include "CpptrajStdio.h"

// CONSTRUCTOR
FileIO_Bzip2::FileIO_Bzip2() {
  //fprintf(stderr,"FileIO_Bzip2 CONSTRUCTOR\n");
  fp_ = NULL;
  infile_ = NULL;
  err_ = BZ_OK;
  isBzread_ = true;
  bzfilename_ = NULL;
  bzmode_ = NULL;
  position_ = 0L;
}

// DESTRUCTOR
FileIO_Bzip2::~FileIO_Bzip2() {
  //fprintf(stderr,"FileIO_Bzip2 DESTRUCTOR\n");
  if (fp_!=NULL || infile_!=NULL) this->Close();
  if (bzfilename_!=NULL) free(bzfilename_);
  if (bzmode_!=NULL) free(bzmode_);
}

// FileIO_Bzip2::BZerror()
/** Return a string corresponding to the current value of err.
  */
const char *FileIO_Bzip2::BZerror() {
  switch (err_) {
    case BZ_OK : return "BZ_OK";
    case BZ_PARAM_ERROR : return "BZ_PARAM_ERROR";
    case BZ_SEQUENCE_ERROR : return "BZ_SEQUENCE_ERROR";
    case BZ_IO_ERROR: return "BZ_IO_ERROR";
    case BZ_UNEXPECTED_EOF : return "BZ_UNEXPECTED_EOF";
    case BZ_DATA_ERROR : return "BZ_DATA_ERROR";
    case BZ_DATA_ERROR_MAGIC : return "BZ_DATA_ERROR_MAGIC";
    case BZ_MEM_ERROR : return "BZ_MEM_ERROR";
    case BZ_STREAM_END : return "BZ_MEM_ERROR";
  }
  return "Unknown Bzip2 error";
}

// FileIO_Bzip2::Open()
/** Open the given file as a bzip2 file. The mode and filename are stored
  * in case rewind is called (Bzip2 routines do not have a rewind so the
  * file must be closed and reopened).
  */
// NOTES from Bzip2 docs:
// BZFILE *BZ2_bzReadOpen( int *bzerror, FILE *f, int verbosity, int small,
//                         void *unused, int nUnused );
//   If small is 1, the library will try to decompress using less memory, at 
//   the expense of speed. BZ2_bzRead will decompress the nUnused bytes starting 
//   at unused, before starting to read from the file f. At most BZ_MAX_UNUSED 
//   bytes may be supplied like this. If this facility is not required, you 
//   should pass NULL and 0 for unused and nUnused respectively.
// BZFILE *BZ2_bzWriteOpen( int *bzerror, FILE *f, int blockSize100k, 
//                          int verbosity, int workFactor );
//   Parameter blockSize100k specifies the block size to be used for compression.
//   It should be a value between 1 and 9 inclusive, and the actual block size 
//   used is 100000 x this figure. 9 gives the best compression but takes most 
//   memory. Parameter verbosity should be set to a number between 0 and 4 
//   inclusive. 0 is silent. Parameter workFactor controls how the compression 
//   phase behaves when presented with worst case, highly repetitive, input data.
//   If compression runs into difficulties caused by repetitive data, the library 
//   switches from the standard sorting algorithm to a fallback algorithm. The 
//   fallback is slower than the standard algorithm by perhaps a factor of three, 
//   but always behaves reasonably, no matter how bad the input. Lower values of 
//   workFactor reduce the amount of effort the standard algorithm will expend 
//   before resorting to the fallback. You should set this parameter carefully; 
//   too low, and many inputs will be handled by the fallback algorithm and so 
//   compress rather slowly, too high, and your average-to-worst case compression
//   times can become very large. The default value of 30 gives reasonable 
//   behaviour over a wide range of circumstances.
int FileIO_Bzip2::Open(const char *filename, const char *mode) {
  // Store filename and mode - reallocate in case of reopen
  if (bzfilename_!=filename) {
    bzfilename_ = (char*) realloc(bzfilename_, (strlen(filename)+1) * sizeof(char));
    strcpy(bzfilename_, filename);
  }
  if (bzmode_!=mode) {
    bzmode_     = (char*) realloc(bzmode_,     (strlen(mode)+1    ) * sizeof(char));
    strcpy(bzmode_, mode);
  }

  // DEBUG
  //mprintf("DEBUG: FileIO_Bzip2::Open(%s,%s)\n",filename,mode);

  fp_ = fopen(filename, mode);
  if (fp_==NULL) {
    mprintf("Error: FileIO_Bzip2::Open: Could not open %s with mode %s\n",filename,mode);
    return 1;
  }

  switch ( mode[0] ) {
    case 'r' : 
      //mprintf("DEBUG: Calling bzReadOpen\n");
      infile_ = BZ2_bzReadOpen( &err_, fp_, 1, 0, NULL, 0); 
      isBzread_ = true;  
      break;
    case 'w' : 
      //mprintf("DEBUG: Calling bzWriteOpen\n");
      infile_ = BZ2_bzWriteOpen( &err_, fp_, 9, 0, 30);
      isBzread_ = false; 
      break;
    case 'a' : 
      mprintf("Error: FileIO_Bzip2::Open: Append not supported for Bzip2.\n");
      return 1; // No append for Bzip2
    default: return 1; 
  }

  if (err_ != BZ_OK) {
    mprintf("Error: FileIO_Bzip2::Open: Could not BZOPEN %s with mode %s\n",filename,mode);
    return 1;
  }

  if (infile_==NULL) return 1;
  //mprintf("DEBUG: BZIP2 Opened %s with mode %s\n",filename,mode);
  position_ = 0L;
  return 0;
}

// FileIO_Bzip2::Close()
int FileIO_Bzip2::Close() {
  if (infile_!=NULL) {
    if (isBzread_) {
      //mprintf("DEBUG: BZ2_bzReadClose\n");
      BZ2_bzReadClose(&err_, infile_);
    } else {
      //mprintf("DEBUG: BZ2_bzWriteClose\n");
      BZ2_bzWriteClose(&err_, infile_, 0, NULL, NULL);
    }
    infile_ = NULL;
  }
  
  if (fp_!=NULL) fclose(fp_);
  fp_ = NULL;
  return 0;
}

// FileIO_Bzip2::Size()
/** Since the uncompressed size of Bzip files is not stored anywhere in the file
  * need to read every possible byte in the file, which can be VERY slow. 
  * NOTE: The input filename is currently IGNORED.
  * NOTE: This can be ridiculously time consuming for large bzip files, so
  *       just return 0. 
  */
//#define BUFINSIZE 10240
off_t FileIO_Bzip2::Size(const char *filename) {
  //off_t fileSize, numread;
  //char bufIn[BUFINSIZE];
//  char Scan;

  if (filename==NULL) return -1L;
  //fileSize=0L;
  return 0L;
/*
  // Check that the file being checked is the currently open file.
  // NOTE: Unnecessary?
  //if (strcmp(filename, bzfilename)!=0) {
  //  mprintf("ERROR: FileIO_Bzip2::Size: Checking file %s, open file is %s!\n",
  //          filename,bzfilename);
  //  return -1L;
  //}

  // Open the file
  if (infile==NULL) {
    mprintf("FileIO_Bzip2::Size: Opening %s\n",filename);
    if (this->Open(filename,"rb")) return -1L;
  }

  // Read all chars in file 10K bytes at a time.
  // NOTE: Use sizeof(char)??
  // NOTE: Larger buffer? Dynamically allocate?
//  bufIn = (char*) malloc(10240 * sizeof(char));
  while ( (numread = (off_t) this->Read(bufIn, 1, BUFINSIZE)) > 0 )
    fileSize += numread;
//  free(bufIn);
//  while ( this->Read(&Scan, 1, 1) > 0 )
//    fileSize = fileSize + 1L;

  // Close file
  this->Close();

  //mprintf("FileIO_Bzip2::Size: Uncompressed size of %s: %lu\n",filename,fileSize);

  return fileSize;
*/
}
//#undef BUFINSIZE

// FileIO_Bzip2::Read()
/** Read num_bytes from bzip2file stream.  
  * \return number of bytes read on success.
  * \return -1 on error.
  */
int FileIO_Bzip2::Read(void *buffer, size_t num_bytes) {
  int numread = BZ2_bzRead(&err_, infile_, buffer, num_bytes);
  // Update position
  position_ += ((off_t) numread);
  if (err_!=BZ_OK && err_!=BZ_STREAM_END) {
    mprinterr("Error: FileIO_Bzip2::Read: BZ2_bzRead error: [%s]\n"
              "Error:                     size=%i expected=%zu\n",
               this->BZerror(), numread, num_bytes);
    return -1;
  }
  return numread;
}

// FileIO_Bzip2::Write()
int FileIO_Bzip2::Write(const void *buffer, size_t num_bytes) {
  // NOTE: The bzip2 library requires the void* cast
  BZ2_bzWrite ( &err_, infile_, (void*)buffer, num_bytes );
  // Update position
  position_ += ((off_t)num_bytes);
  if (err_ == BZ_IO_ERROR) { 
    mprintf( "Error: FileIO_Bzip2::Write: BZ2_bzWrite error\n");
    return 1;
  }
  return 0;
}

// FileIO_Bzip2::Seek
/** Since a true seek is not really possible with bzip2, scan 1 char at 
  * a time until desired position achieved.
  */
// NOTE: Scan in blocks?
int FileIO_Bzip2::Seek(off_t offset) {
  off_t seekTo;
  char Scan;
  // Determine place to seek to
  //switch (origin) {
  //  case SEEK_SET : seekTo = offset; break;
  //  case SEEK_CUR : seekTo = position + offset; break;
  //  case SEEK_END : 
  //    mprintf("Error: FileIO_Bzip2::Seek: Seek to END not supported (%s).\n",bzfilename);
  //  default : return 1;
  //}
  seekTo = offset;

  //mprintf("DEBUG: FileIO_Bzip2::Seek: %s %li -> %li, ",bzfilename,position,seekTo);

  // If place to seek to is earlier than current position need to reopen
  if (seekTo<position_) 
    this->Rewind();

  // Read chars until position achieved
  while (position_ < seekTo) {
    if (this->Read(&Scan,1) < 1) break;
  }

  //mprintf("%li\n",position);

  return 0;
}

// FileIO_Bzip2::Rewind()
/** Close and reopen.
  */
int FileIO_Bzip2::Rewind() {
  if (bzfilename_==NULL || bzmode_==NULL) return 1;
  this->Close();
  this->Open(bzfilename_,bzmode_);
  return 0;
}

// FileIO_Bzip2::Tell()
// NOTE: Tell not possible with bzip2. Use position.
off_t FileIO_Bzip2::Tell() {
  return position_;
}

// FileIO_Bzip2::Gets()
/** Analogous to fgets, reads characters from stream and stores them as a C 
  * string into str until (num-1) characters have been read or either a newline
  * or the End-of-File is reached, whichever comes first.
  * A newline character makes fgets stop reading, but it is considered a valid 
  * character and therefore it is included in the string copied to str.
  * A null character is automatically appended in str after the characters read
  * to signal the end of the C string.
  */
int FileIO_Bzip2::Gets(char *str, int num) {
  int i;
  //mprintf("DEBUG: FileIO_Bzip2::Gets: num=%i\n",num);
  // Try to read num chars. If newline encountered exit
  if (num<=1) return 1;
  i=0;
  while ( this->Read(str+i, 1) > 0 ) {
    i++;
    if (i==num-1) break;
    if (str[i-1]=='\n') break;
  }
  // If nothing read return 1
  if (i==0) return 1;
  // i should be at num or 1 after newline; append NULL char
  str[i] = '\0';
  //mprintf("DEBUG: FileIO_Bzip2::Gets: num=%i i=%i [%s]\n",num,i,str);
  //mprintf( "DEBUG: After FileIO_Bzip2::Gets: position %li\n",position);
  return 0;
}
#endif
