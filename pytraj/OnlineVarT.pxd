# distutils: language = c++
from libcpp.map cimport map

cdef extern from "OnlineVarT.h": 
    ctypedef map[int, double].iterator iterator
    cdef cppclass _Stats "Stats" [Float]:
        _Stats()
        void accumulate(const Float)
        Float mean()
        Float variance()
        Float nData()
    cdef cppclass _StatsMap "StatsMap" [int, double]:
        _StatsMap() 
        void accumulate(map[int, double] a)
        double mean(int i)
        double variance(int i)
        iterator mean_begin() 
        iterator mean_end() 
        iterator variance_begin() 
        iterator variance_end() 
        double nData() const 

cdef class StatsMap:
    cdef _StatsMap* thisptr

cdef class Stats:
    cdef _Stats* thisptr
