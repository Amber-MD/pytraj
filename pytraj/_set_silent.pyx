# distutils: language = c++

cdef extern from "CpptrajStdio.h":
    void SupressErrorMsg(bint)

def set_error_silent(turnoff=True):
    SupressErrorMsg(turnoff)
