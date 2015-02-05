// FileIO_Std: Standard C file operations
#include "FileIO_Std.h" // FileIO.h, cstdio

// CONSTRUCTOR
FileIO_Std::FileIO_Std() {
  fp_ = NULL;
  isStdout_=false;
}

// DESTRUCTOR
FileIO_Std::~FileIO_Std() {
  if (fp_!=NULL) this->Close();
}

// FileIO_Std::Open()
/** Open file using standard C routines. If mode is WRITE and no
  * filename given default to stdout.
  */
int FileIO_Std::Open(const char *filename, const char *mode) {
  if (filename==NULL) {
    if (mode[0]=='w') 
      fp_=stdout;
    else
      return 1;
    isStdout_=true;
  } else
    fp_ = fopen(filename, mode);
  if (fp_==NULL) return 1;
  return 0;
}

// FileIO_Std::Close()
/** Close stream if not stdout
  */
int FileIO_Std::Close() {
  if (fp_!=NULL && !isStdout_) fclose(fp_);
  fp_=NULL;
  return 0;
}

// FileIO_Std::Read()
int FileIO_Std::Read(void *buffer, size_t num_bytes) {
  size_t numread = fread(buffer, 1, num_bytes, fp_);
  if (ferror(fp_)) {
    perror("Error during FileIO_Std::Read");
    return -1;
  }
  return (int) numread;
}

// FileIO_Std::Write()
int FileIO_Std::Write(const void *buffer, size_t num_bytes) {
  size_t numwrite = fwrite(buffer, 1, num_bytes, fp_);
  // NOTE: Check for errors here.
  if (numwrite != num_bytes) return 1;
  return 0;
}

// FileIO_Std::Seek()
// NOTE: Use fseeko for better compatibility with large files.
int FileIO_Std::Seek(off_t offset) {
  // DEBUG
  //printf("Calling standard seek(%i): %li\n",origin,offset);
#ifdef _MSC_VER
	return _fseeki64(fp_,offset,SEEK_SET);
#else
  return fseeko(fp_, offset, SEEK_SET);
#endif
}

// FileIO_Std::Rewind()
int FileIO_Std::Rewind() {
  rewind(fp_);
  return 0;
}

// FileIO_Std::Tell()
off_t FileIO_Std::Tell() {
#ifdef _MSC_VER
	return _ftelli64(fp_);
#else
  return ftello(fp_);
#endif
}

// FileIO_Std::Gets()
int FileIO_Std::Gets(char *str, int num) {
  if ( fgets(str,num,fp_) == NULL ) {
    //fprintf(stdout,"DEBUG: FileIO_Std::Gets returned NULL (%s) %i\n",str,num);
    return 1;
  } else
    return 0;
}

