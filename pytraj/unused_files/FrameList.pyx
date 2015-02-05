# distutils: language = c++


cdef class FrameList:
    """Fake class
    TODO : cpptraj (> v15.22b) no longer support this class
    """
    def __cinit__(self):
        self.thisptr = new _FrameList()

    def __dealloc__(self):
        del self.thisptr
