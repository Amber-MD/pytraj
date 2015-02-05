// FileIO_Mpi
#include <cstdio>
#include <cstdlib>
#include "FileIO_Mpi.h"

// CONSTRUCTOR
FileIO_Mpi::FileIO_Mpi() {
  pfile_ = (parallelType) malloc( sizeof(parallelType));
}

// DESTRUCTOR
FileIO_Mpi::~FileIO_Mpi() {
  free(pfile_);
}

// FileIO_Mpi::Open()
int FileIO_Mpi::Open(const char *filename, const char *mode) {
  int err=0;

  switch( mode[0] ) {
    case 'r' : err=parallel_openFile_read(pfile_, filename); break;
    case 'w' : err=parallel_open_file_write(pfile_, filename); break;
    case 'a' : err=1; break; // NOTE: No MPI append for now
    default  : err=1; break;
  }

  return err;
}

// FileIO_Mpi::Close()
int FileIO_Mpi::Close() {
  parallel_closeFile(pfile_);
  return 0;
}

// FileIO_Mpi::Read()
int FileIO_Mpi::Read(void *buffer, size_t num_bytes) {
  return (parallel_fread(pfile_, buffer, num_bytes));
}

// FileIO_Mpi::Write()
int FileIO_Mpi::Write(const void *buffer, size_t num_bytes) {
  if ( parallel_fwrite(pfile_, buffer, num_bytes) ) return 1;
  // NOTE: Check for errors here.
  return 0;
}

int FileIO_Mpi::Flush() {
  return parallel_flush( pfile_ );
}

// FileIO_Mpi::Seek()
int FileIO_Mpi::Seek(off_t offset) {

  if ( parallel_fseek(pfile_, offset, SEEK_SET) ) return 1;

  return 0;
}

// FileIO_Mpi::Rewind()
int FileIO_Mpi::Rewind() {
  if ( parallel_fseek(pfile_, 0L, SEEK_SET) ) return 1;
  return 0;
}

// FileIO_Mpi::Tell()
off_t FileIO_Mpi::Tell() {
  return ( parallel_position(pfile_) );
}

// FileIO_Mpi::Gets()
int FileIO_Mpi::Gets(char *str, int num) {

  if ( parallel_fgets(pfile_,str,num) == NULL ) return 1;
  return 0;
}

// FileIO_Mpi::SetSize()
/** Set size of mpi file, required when splitting up writes.
  */
int FileIO_Mpi::SetSize(long int offset) {

  if ( parallel_setSize(pfile_, offset) ) return 1;
  return 0;
}

