#ifndef INC_BUFFEREDFRAME_H
#define INC_BUFFEREDFRAME_H
#include "CpptrajFile.h"
/// Used to buffer text files that will be read in chunks
class BufferedFrame : public CpptrajFile {
  public:
    BufferedFrame();
    ~BufferedFrame();

    size_t SetupFrameBuffer(int, int, int);
    size_t SetupFrameBuffer(int, int, int, size_t, int);
    size_t ResizeBuffer(int);
    int SeekToFrame(size_t);
    /// Attempt to read frameSize_ bytes.
    int AttemptReadFrame();
    /// Read frameSize_ bytes.
    bool ReadFrame();
    int WriteFrame();
    void GetDoubleAtPosition(double&,size_t,size_t);
    void BufferBegin();
    void BufferBeginAt(size_t);
    void AdvanceBuffer(size_t);
    void BufferToDouble(double*,int);
    void DoubleToBuffer(const double*,int, const char*);

    size_t FrameSize()   const { return frameSize_; }
    const char* Buffer() const { return buffer_;    }
  private:
    size_t CalcFrameSize(int) const;

    char* buffer_;         ///< Character buffer.
    char* bufferPosition_; ///< Position in buffer.
    size_t frameSize_;     ///< Total size of frame to read.
    size_t offset_;        ///< User specified offset, used in seeking.
    int Ncols_;            ///< Number of columns, use to convert array to buffer.
    size_t eltWidth_;      ///< Width of each element in the frame.
};
#endif
