# distutils: language = c++
from libcpp.vector cimport vector


cdef extern from "ClusterSieve.h": 
    ctypedef vector[int] SievedFrames
    # ClusterSieve.h
    ctypedef enum SieveType "ClusterSieve::SieveType":
        NONE "ClusterSieve::NONE"
        REGULAR "ClusterSieve::REGULAR"
        RANDOM "ClusterSieve::RANDOM"
    cdef cppclass _ClusterSieve "ClusterSieve":
        _ClusterSieve() 
        int SetSieve(int, size_t, int)
        int SetSieve(int, const vector[bint]&)
        SievedFrames _Frames() const 
        inline int _FrameToIdx(int frame) const 
        inline size_t Max_Frames() const 
        inline int Sieve() const 
        inline SieveType Type() const 
