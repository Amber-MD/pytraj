# distutils: language = c++
from ClusterList cimport *


cdef extern from "Cluster_HierAgglo.h": 
    # Cluster_HierAgglo.h
    ctypedef enum LINKAGETYPE "Cluster_HierAgglo::LINKAGETYPE":
        SINGLELINK "Cluster_HierAgglo::SINGLELINK"
        AVERAGELINK "Cluster_HierAgglo::AVERAGELINK"
        COMPLETELINK "Cluster_HierAgglo::COMPLETELINK"
    cdef cppclass _Cluster_HierAgglo "Cluster_HierAgglo":
        _Cluster_HierAgglo() 
        void Help() 
        int SetupCluster(_ArgList&)
        void ClusteringInfo() 
        int Cluster() 
        void AddSieved_Frames() 
        void ClusterResults(_CpptrajFile&) const 
