#ifndef INC_CLUSTER_READINFO_H
#define INC_CLUSTER_READINFO_H
#include "ClusterList.h"
/// Read in previous cluster run
class Cluster_ReadInfo : public ClusterList {
  public:
    Cluster_ReadInfo();
    static void Help();
    int SetupCluster(ArgList&);
    void ClusteringInfo();
    int Cluster();
    void AddSievedFrames() { }
    void ClusterResults(CpptrajFile&) const;
  private:
    std::string filename_;
    std::string algorithm_;
};
#endif
