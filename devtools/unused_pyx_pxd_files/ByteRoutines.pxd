# distutil: language = c++

cdef extern from "ByteRoutines.h":
    ctypedef union byte8:
        unsigned char c[8]
        int i[2]
        double d

    void endian_swap(void*, long)
    void endian_swap8(void*, long)
