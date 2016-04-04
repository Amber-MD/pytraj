# distutils: language = c++

from xxx cimport *

cdef class ArgList:
    cdef _yyy  *thisptr

    def __cinit__(self):
        self.thisptr = new _yyy()

    def __dealloc__(self):
        del self.thisptr

