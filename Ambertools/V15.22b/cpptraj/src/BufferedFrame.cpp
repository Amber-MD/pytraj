#include <cstdio>  // sprintf
#include <cstdlib> // atof
#include "BufferedFrame.h"
#include "CpptrajStdio.h"

BufferedFrame::BufferedFrame() :
  buffer_(0),
  bufferPosition_(0),
  frameSize_(0),
  offset_(0),
  Ncols_(0),
  eltWidth_(0)
{}

BufferedFrame::~BufferedFrame() {
  if (buffer_!=0) delete[] buffer_;
}

// BufferedFrame::SetupFrameBuffer()
size_t BufferedFrame::SetupFrameBuffer(int Nelts, int eltWidthIn, int eltsPerLine)
{
  return SetupFrameBuffer(Nelts, eltWidthIn, eltsPerLine, 0, 0);
}

/** Prepare the buffer to receive organized chunks of text, i.e. 
  * organized in some regular fashion (e.g. an Amber Traj, which
  * is 10 cols of 8.3 precision floating point numbers etc).
  * \param Nelts Total expected number of elements to read.
  * \param eltWidth Width in chars of each element.
  * \param eltsPerLine Number of elements per line (columns).
  * \param additionalBytes Any additional bytes in a frame.
  * \param offsetIn Offset (not part of frame) to be used in seeking.
  * \return Size of set-up frame.
  */
size_t BufferedFrame::SetupFrameBuffer(int Nelts, int eltWidthIn, int eltsPerLine, 
                                      size_t additionalBytes, int offsetIn) 
{
  Ncols_ = eltsPerLine;
  eltWidth_ = (size_t)eltWidthIn;
  offset_ = (size_t) offsetIn;
  frameSize_ = CalcFrameSize( Nelts ) + additionalBytes;
  if (buffer_!=0) delete[] buffer_;
  if (frameSize_ < 1) 
    buffer_ = 0;
  else {
    buffer_ = new char[ frameSize_ ];
    std::fill(buffer_, buffer_ + frameSize_, 0);
  }
  bufferPosition_ = buffer_;
  return frameSize_;
}

/** Based on the current values of Ncols and eltWidth, calculate size
  * in bytes necessary to buffer Nelts.
  */
size_t BufferedFrame::CalcFrameSize( int Nelts ) const {
  int frame_lines = Nelts / Ncols_;
  if ((Nelts % Ncols_) > 0) 
    ++frame_lines;
  bool readingFile = (Access() == CpptrajFile::READ);
  // If Reading and DOS, CR present for each newline
  if (readingFile && IsDos()) frame_lines *= 2;
  // Calculate total frame size
  size_t fsize = (((size_t)Nelts * eltWidth_) + frame_lines);
  // If writing, add +1 for null 
  if (!readingFile)
    ++fsize;
  return fsize;
}

/** Increase size of buffer by delta elements, keeping contents intact. */
size_t BufferedFrame::ResizeBuffer(int delta) {
  if (delta == 0) return frameSize_;
  if (delta < 0) {
    mprinterr("Internal Error: ResizeBuffer: Negative value given.\n");
    return frameSize_;
  }
  size_t newsize = frameSize_ + CalcFrameSize( delta );
  char* newbuffer = new char[ newsize ];
  std::copy(buffer_, buffer_+frameSize_, newbuffer);
  std::fill(newbuffer+frameSize_, newbuffer+newsize, 0);
  delete[] buffer_;
  buffer_ = newbuffer;
  bufferPosition_ = buffer_;
  frameSize_ = newsize;
  return frameSize_;
}

int BufferedFrame::SeekToFrame(size_t set) {
  return Seek( (off_t)((set * frameSize_) + offset_) );
}

/** Attempt to read in the next frameSize_ bytes.
  * \return the actual number of bytes read.
  */
int BufferedFrame::AttemptReadFrame() {
  return Read( buffer_, frameSize_ );
}

/** Read the next frameSize_ bytes.
  * \return true if the number of bytes read was less than frameSize_.
  * \return false if read was successful.
  */
bool BufferedFrame::ReadFrame() {
  return ( Read( buffer_, frameSize_ ) != (int)frameSize_ );
}

int BufferedFrame::WriteFrame() {
  return Write( buffer_, (size_t)(bufferPosition_ - buffer_) );
}

void BufferedFrame::GetDoubleAtPosition(double& val, size_t start, size_t end) {
  char savechar = buffer_[end];
  buffer_[end] = '\0';
  val = atof(buffer_ + start);
  buffer_[end] = savechar;
}

void BufferedFrame::BufferBegin() {
  bufferPosition_ = buffer_;
}

void BufferedFrame::BufferBeginAt(size_t pos) {
  bufferPosition_ = buffer_ + pos;
}

void BufferedFrame::AdvanceBuffer(size_t offset) {
  bufferPosition_ += offset;
}

/** Convert text in buffer containing numerical elements with format 
  * X0Y0Z0X1Y1Z1...XNYNZN to the given double array. The width of each 
  * element should be what SetupFrameBuffer was called with, and the 
  * number of elements to read should not be greater than Nelts.
  * Newlines are skipped. Output array should be as big as Nout. 
  * Update bufferPosition after read.
  */
void BufferedFrame::BufferToDouble(double* Xout, int Nout) {
  for (int element = 0; element < Nout; ++element) {
    // Advance past newlines / CR (dos)
    while (*bufferPosition_=='\n' || *bufferPosition_=='\r')
      ++bufferPosition_;
    if (*bufferPosition_ == '*') {
      mprinterr("Error: '*' encountered (atom %i", (element / 3) + 1);
      int problem_xyz = element % 3;
      if (problem_xyz == 0)      mprinterr(" X");
      else if (problem_xyz == 1) mprinterr(" Y");
      else                       mprinterr(" Z");
      mprinterr("). This indicates coordinate overflow.\n");
    }
    char *ptrend = bufferPosition_ + eltWidth_;
    char lastchar = *ptrend;
    *ptrend = '\0';
    Xout[element] = atof(bufferPosition_);
    *ptrend = lastchar;
    bufferPosition_ = ptrend;
  }
}

/** Convert given double array to ordered text in buffer. The number of
  * elements in the given array should be what SetupFrameBuffer was
  * called with. Update bufferPosition after write. 
  */
void BufferedFrame::DoubleToBuffer(const double* Xin, int Nin, const char* format)
{
  int col = 0;
  for (int element = 0; element < Nin; ++element) {
    sprintf(bufferPosition_, format, Xin[element]);
    bufferPosition_ += eltWidth_;
    ++col;
    if ( col == Ncols_ ) {
      sprintf(bufferPosition_,"\n");
      ++bufferPosition_;
      col = 0;
    }
  }
  // If the coord record didnt end on a newline, print one
  if ( col != 0 ) {
    sprintf(bufferPosition_,"\n");
    ++bufferPosition_;
  }
}
