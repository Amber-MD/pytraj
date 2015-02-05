#ifndef INC_CLUSTERLIST_H
#define INC_CLUSTERLIST_H
#include <list>
#include "ArgList.h"
#include "ClusterNode.h" 
// Class: ClusterList
/** This base class holds all the individual clusters, as well as routines 
  * that can be used to obtain information on clusters after clustering.
  */
class ClusterList {
  public:
    enum DistModeType   { USE_FRAMES = 0, USE_FILE  };
    enum DistMetricType { RMS = 0, DME, SRMSD, DATA };
    static const char* MetricString( DistMetricType );
    ClusterList();
    virtual ~ClusterList();
    int Nclusters()                  const { return (int)clusters_.size(); }

    void SetDebug(int);
    void Renumber(bool);
    void Summary(std::string const&,int);
    void Summary_Part(std::string const&,int,std::vector<int> const&);
    void PrintClustersToFile(std::string const&,int);
    void PrintClusters();

    int CalcFrameDistances(std::string const&, ClusterDist::DsArray const&, DistModeType, 
                           DistMetricType, bool, bool, std::string const&, int, int);
    // Inherited by individual clustering methods
    virtual int SetupCluster(ArgList&) = 0;
    virtual void ClusteringInfo() = 0;
    virtual int Cluster() = 0;

    // Const Iterator over clusters
    typedef std::list<ClusterNode>::const_iterator cluster_iterator;
    const cluster_iterator begincluster() const { return clusters_.begin(); }
    const cluster_iterator endcluster()   const { return clusters_.end();   }
    /// Remove clusters with no members.
    void RemoveEmptyClusters();
    /// Calculate distances between each cluster
    void CalcClusterDistances();
    /// Calculate cluster silhouettes
    void CalcSilhouette(std::string const&) const;

    void DrawGraph(bool,DataSet*,double,int) const;
  protected:
    virtual void AddSievedFrames() = 0;
    virtual void ClusterResults(CpptrajFile&) const = 0;

    void AddSievedFramesByCentroid();
    /// Iterator over clusters
    typedef std::list<ClusterNode>::iterator cluster_it;
    int debug_;
    /// Store individual cluster info; frame numbers, centroid, etc.
    std::list<ClusterNode> clusters_;
    /// Distances between each frame.
    ClusterMatrix FrameDistances_;
    /// Distances between each cluster.
    ClusterMatrix ClusterDistances_;
    /// Used to calculate distances between frames and/or centroids.
    ClusterDist* Cdist_;
    /// Add specified frames to a new cluster.
    int AddCluster(ClusterDist::Cframes const&);
    /// Calculate the Davies-Bouldin index of clusters.
    double ComputeDBI(CpptrajFile&);
    /// Calculate pseudo-F statistic.
    double ComputePseudoF(CpptrajFile&);
  private:
    static const char* XMGRACE_COLOR[];
};
#endif
