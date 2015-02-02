#ifndef INC_CLUSTERNODE_H
#define INC_CLUSTERNODE_H
#include "ClusterDist.h" 
/// Hold information for a cluster in a ClusterList
class ClusterNode {
  public:
    ClusterNode();
    ~ClusterNode();
    ClusterNode(ClusterDist*,ClusterDist::Cframes const&, int);
    ClusterNode(const ClusterNode&);
    ClusterNode& operator=(const ClusterNode&);
    /// Used to sort clusters by # of frames in cluster
    inline bool operator<(const ClusterNode&) const;
    /// Merge frames from another cluster to this cluster
    inline void MergeFrames(ClusterNode const&);
    /// Determine which frame in the cluster is the best representative.
    int FindBestRepFrame(ClusterMatrix const&);
    /// Calculate eccentricity for frames in this cluster.
    void CalcEccentricity(ClusterMatrix const&);
    /// Calculate centroid of members of this cluster.
    void CalculateCentroid(ClusterDist* Cdist) {
      // FIXME: Could potentially get rid of this branch.
      if (centroid_ == 0)
        centroid_ = Cdist->NewCentroid( frameList_ );
      else
        Cdist->CalculateCentroid( centroid_, frameList_ );
    }
    /// Calculate average distance of all members to centroid
    double CalcAvgToCentroid( ClusterDist*);
    // Iterator over frame numbers
    typedef ClusterDist::Cframes::const_iterator frame_iterator;
    frame_iterator beginframe() const { return frameList_.begin(); }
    frame_iterator endframe()   const { return frameList_.end();   }
    int ClusterFrame(int idx)         const { return frameList_[idx];    } 
    // Return internal variables
    inline double AvgDist()      const { return avgClusterDist_;        }
    inline double Eccentricity() const { return eccentricity_;          }
    inline int Num()             const { return num_;                   }
    inline int Nframes()         const { return (int)frameList_.size(); }
    inline int BestRepFrame()   const  { return repFrame_;              }
    inline Centroid* Cent()            { return centroid_;              }
    // Set internal variables 
    void SetAvgDist(double avg)        { avgClusterDist_ = avg;         }
    void AddFrameToCluster(int fnum)   { frameList_.push_back( fnum );  }
    void SetNum(int numIn)             { num_ = numIn;                  }
    void SortFrameList();
    bool HasFrame(int);
    void RemoveFrameFromCluster(int);
    void RemoveFrameUpdateCentroid(ClusterDist*, int);
    void AddFrameUpdateCentroid(ClusterDist*, int);
  private:
    double avgClusterDist_;           ///< Avg distance of this cluster to all other clusters.
    double eccentricity_;             ///< Maximum distance between any 2 frames.
    int num_;                         ///< Cluster number.
    int repFrame_;                    ///< Frame number with lowest dist. to all other frames.
    ClusterDist::Cframes frameList_;  ///< List of frames belonging to this cluster.
    Centroid* centroid_;              ///< Centroid of all frames in this cluster. 
};
// ----- INLINE FUNCTIONS ------------------------------------------------------
/** Use > since we give higher priority to larger clusters. */
bool ClusterNode::operator<(const ClusterNode& rhs) const {
  return ( frameList_.size() > rhs.frameList_.size() );
}
/** Frames from rhs go to this cluster. */
void ClusterNode::MergeFrames( ClusterNode const& rhs) {
  frameList_.insert(frameList_.end(), rhs.frameList_.begin(), rhs.frameList_.end());
}
#endif
