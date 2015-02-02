# distutil: language = c++

cdef import from "CpptrajStdio.h":
    void mflush()
    void mprintf(const char*, ...)
    void mprinterr(const char*, ...)
    void rprintf(const char*, ...)
    void rprinterr(const char*, ...)
    void SetWorldSilent(bint)
