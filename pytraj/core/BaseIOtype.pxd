# distutil: language = c++

ctypedef _BaseIOtype* (*AllocatorType)()
ctypedef void (*HelpType)()

cdef extern from "BaseIOtype.h":
    #ctypedef _BaseIOtype* (*AllocatorType)()
    #ctypedef void (*HelpType)()
    cdef cppclass _BaseIOtype "BaseIOtype":
        pass

cdef class BaseIOtype:
    cdef _BaseIOtype* baseptr0
