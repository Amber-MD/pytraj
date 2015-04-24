# distutil: language = c++

cdef extern from "FileIO.h":
    cdef cppclass _FileIO "FileIO":
        pass
