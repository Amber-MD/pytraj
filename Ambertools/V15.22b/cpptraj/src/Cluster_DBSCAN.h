#ifndef INC_CLUSTER_DBSCAN_H
#define INC_CLUSTER_DBSCAN_H
#include "ClusterList.h"
class Cluster_DBSCAN : public ClusterList {
  public:
    Cluster_DBSCAN();
    static void Help();
    int SetupCluster(ArgList&);
    void ClusteringInfo();
    int Cluster();
    void AddSievedFrames();
    void ClusterResults(CpptrajFile&) const;
  private:
    static char UNASSIGNED;
    static char NOISE;
    static char INCLUSTER;
    int minPoints_;            ///< Min # of points needed to make a cluster.
    double epsilon_;           ///< Distance criterion for cluster formation.
    Range kdist_;
    std::string k_prefix_;     ///< Kdist output file prefix.
    bool sieveToCentroid_;     ///< If true sieve only based on closeness to centroid.
    std::vector<char> Status_; ///< Hold current status of each frame.
    void RegionQuery(std::vector<int>&, std::vector<int> const&, int);
    void ComputeKdist( int, std::vector<int> const& ) const ;
    void ComputeKdistMap( Range const&, std::vector<int> const& ) const ;
};
#endif
