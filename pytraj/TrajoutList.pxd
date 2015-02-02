# distutils: language = c++
from libcpp.vector cimport vector
from libcpp.string cimport string
from pytraj.Trajout cimport *
from pytraj.TopologyList cimport *


cdef extern from "TrajoutList.h": 
    cdef cppclass _TrajoutList "TrajoutList":
        _TrajoutList() 
        #~_TrajoutList() 
        void Clear() 
        void SetDebug(int)
        int AddEnsembleTrajout(const _ArgList&, const _TopologyList&, int)
        int AddTrajout(const _ArgList&, const _TopologyList&)
        int WriteTrajout(int, _Topology *, _Frame *)
        void CloseTrajout() 
        void List() const 
        bint Empty() const 
        ArgIt argbegin() const 
        ArgIt argend() const 

cdef class TrajoutList:
    cdef _TrajoutList* thisptr

