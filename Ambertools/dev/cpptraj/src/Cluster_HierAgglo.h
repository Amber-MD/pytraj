#ifndef INC_CLUSTER_HIERAGGLO_H
#define INC_CLUSTER_HIERAGGLO_H
#include "ClusterList.h"
class Cluster_HierAgglo : public ClusterList {
  public:
    /// Type of distance calculation between clusters.
    enum LINKAGETYPE  { SINGLELINK = 0, AVERAGELINK, COMPLETELINK };
    Cluster_HierAgglo();
    static void Help();
    int SetupCluster(ArgList&);
    void ClusteringInfo();
    int Cluster();
    void AddSievedFrames() { AddSievedFramesByCentroid(); }
    void ClusterResults(CpptrajFile&) const;
  private:
    int nclusters_;       ///< Target # of clusters.
    double epsilon_;      ///< Once the min distance between clusters is > epsilon, stop.
    LINKAGETYPE linkage_; ///< Cluster Linkage type.
    CpptrajFile eps_v_n_; ///< Write epsilon vs # clusters.

    void InitializeClusterDistances();
    int MergeClosest();
    void calcMinDist(cluster_it&);
    void calcMaxDist(cluster_it&);
    void calcAvgDist(cluster_it&);
};
#endif
