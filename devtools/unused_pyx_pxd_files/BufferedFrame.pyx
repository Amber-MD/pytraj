# distutils: language = c++


cdef class BufferedFrame:
    def __cinit__(self):
        self.thisptr = new _BufferedFrame()

    def __dealloc__(self):
        del self.thisptr

    def setup_frame_buffer(self, int Nelts, int eltWidthIn, int eltsPerLine, 
                           size_t additionalBytes=0, int offset=0):
        self.thisptr.SetupFrameBuffer(Nelts, eltWidthIn, eltsPerLine, additionalBytes, offset)

    def resize_buffer(self, int delta):
        return self.thisptr.ResizeBuffer(delta)

    def seek_to_frame(self,size_t idx):
        return self.thisptr.SeekToFrame(idx)

    def attempt_read_frame(self):
        return self.thisptr.AttemptReadFrame()

    def read_frame(self):
        return self.thisptr.ReadFrame()

    def write_frame(self):
        return self.thisptr.WriteFrame()

    def get_double_at_position(self, double val, size_t start, size_t end):
        self.thisptr.GetDoubleAtPosition(val, start, end)

    def buffer_begin(self):
        self.thisptr.BufferBegin()

    def buffer_begin_at(self,size_t idx):
        self.thisptr.BufferBeginAt(idx)

    def advance_buffer(self,size_t idx):
        self.thisptr.AdvanceBuffer(idx)

    def buffer_to_double(self, double[:] Xout, int Nout):
        self.thisptr.BufferToDouble(&Xout[0], Nout)

    def double_to_buffer(self, double[:] Xin, int Nin, char * _format):
        self.thisptr.DoubleToBuffer(&Xin[0], Nin, _format)

    def frame_size(self):
        return self.thisptr.FrameSize()

    def buffer(self):
        return self.thisptr.Buffer()
