# distutils: language = c++
from ClusterList cimport *


cdef extern from "Cluster_DBSCAN.h": 
    cdef cppclass _Cluster_DBSCAN "Cluster_DBSCAN" (_ClusterList):
        _Cluster_DBSCAN() 
        void Help() 
        int SetupCluster(_ArgList&)
        void ClusteringInfo() 
        int Cluster() 
        void AddSievedFrames() 
        void ClusterResults(_CpptrajFile&) const 


cdef class Cluster_DBSCAN (ClusterList):
    cdef _Cluster_DBSCAN* ptr
