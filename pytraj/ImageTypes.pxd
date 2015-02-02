# distutils: language = c++
from libcpp.vector cimport vector

cdef extern from "ImageTypes.h" namespace "Image":
    ctypedef vector[int] PairType
    ctypedef enum Mode:
        BYMOL
        BYRES
        BYATOM
    inline const char* ModeString(Mode m)
