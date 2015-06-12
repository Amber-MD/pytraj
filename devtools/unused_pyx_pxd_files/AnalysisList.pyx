# distutils: language = c++

cdef class AnalysisList:
    def __cinit__(self):
        self.thisptr = new _AnalysisList()

    def __dealloc__(self):
        del self.thisptr
