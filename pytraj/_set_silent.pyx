# distutils: language = c++

cdef extern from "CpptrajStdio.h":
    void SupressErrorMsg(bint)
    void SetWorldSilent(bint)

def set_error_silent(turnoff=True):
    SupressErrorMsg(turnoff)

def set_world_silent(turnoff=True):
    SetWorldSilent(turnoff)
