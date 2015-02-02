# distutils: language = c++

cdef class DataIO_Mdout:
    cdef _DataIO_Mdout  *thisptr

    def __cinit__(self):
        self.thisptr = new _DataIO_Mdout()

    def __dealloc__(self):
        del self.thisptr

