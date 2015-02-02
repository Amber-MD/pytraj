# distutils: language = c++

cdef class BufferedLine:
    #cdef _BufferedLine  *thisptr

    def __cinit__(self):
        self.thisptr = new _BufferedLine()

    def __dealloc__(self):
        del self.thisptr

    def Line(self):
        return self.thisptr.Line()
