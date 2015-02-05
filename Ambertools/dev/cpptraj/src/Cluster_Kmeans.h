#ifndef INC_CLUSTER_KMEANS_H
#define INC_CLUSTER_KMEANS_H
#include "ClusterList.h"
#include "Random.h"
class Cluster_Kmeans : public ClusterList {
  public:
    Cluster_Kmeans();
    static void Help();
    int SetupCluster(ArgList&);
    void ClusteringInfo();
    int Cluster();
    void AddSievedFrames() { AddSievedFramesByCentroid(); }
    void ClusterResults(CpptrajFile&) const;
  private:
    typedef std::vector<int> Iarray;
    enum KmeansModeType { SEQUENTIAL, RANDOM };

    int FindKmeansSeeds();
    void ShufflePoints(Iarray&);

    Random_Number RN_;
    int nclusters_; ///< Target number of clusters.
    int kseed_;
    int maxIt_;
    Iarray SeedIndices_;
    Iarray FramesToCluster_;
    KmeansModeType mode_;
    bool clusterToClusterCentroid_;
};
#endif
