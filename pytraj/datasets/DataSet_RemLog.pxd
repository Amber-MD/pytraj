# distutils: language = c++
from libcpp.vector cimport vector
from libcpp.map cimport map
from .DataSet cimport *

ctypedef map[double, int] TmapType
cdef extern from "DataSet_RemLog.h": 
    cdef cppclass _DataSet_RemLog "DataSet_RemLog":
        _DataSet_RemLog() 
        _DataSet * Alloc() 
        void AllocateReplicas(int)
        void AddRepFrame(int rep, const _ReplicaFrame& frm)
        const _ReplicaFrame& RepFrame(int exch, int rep) const 
        int NumExchange() const 
        bint ValidEnsemble() const 
        void TrimLastExchange() 
        size_t Size() const 
        int Sync() 
        void Info() const 
        void Add(size_t, const void *)


    cdef cppclass _ReplicaFrame "DataSet_RemLog::ReplicaFrame":
        _Replica_Frame() 
        int SetTremdFrame(const char *, const TmapType&)
        int SetHremdFrame(const char *, const vector[int]&)
        int ReplicaIdx() const 
        int PartnerIdx() const 
        int CoordsIdx() const 
        bint Success() const 
        double Temp0() const 
        double PE_X1() const 
        double PE_X2() const 


cdef class DataSet_RemLog:
    cdef _DataSet_RemLog* thisptr

cdef class ReplicaFrame:
    cdef _ReplicaFrame* thisptr

