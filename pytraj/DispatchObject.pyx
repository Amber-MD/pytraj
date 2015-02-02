# distutils: language = c++


cdef class DispatchObject:
    def __cinit__(self):
        self.thisptr = new _DispatchObject()

    def __dealloc__(self):
        del self.thisptr
