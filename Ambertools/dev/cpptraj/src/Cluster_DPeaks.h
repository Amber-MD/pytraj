#ifndef INC_CLUSTER_DPEAKS_H
#define INC_CLUSTER_DPEAKS_H
#include "ClusterList.h"
class Cluster_DPeaks : public ClusterList {
  public:
    Cluster_DPeaks();
    static void Help();
    int SetupCluster(ArgList&);
    void ClusteringInfo();
    int Cluster();
    void AddSievedFrames();
    void ClusterResults(CpptrajFile&) const;
  private:
    void AssignClusterNum(int, int&);

    std::string dpeaks_;
    std::string rafile_;
    std::string radelta_;
    double epsilon_;
    int avg_factor_;
    bool calc_noise_;
    class Cpoint {
      public:
        Cpoint() : dist_(-1.0), density_(0), fnum_(-1), nidx_(-1), oidx_(-1), cnum_(-1) {}
        //Cpoint(int f, int o) : dist_(-1.0), density_(0), fnum_(f), nidx_(-1), oidx_(o), cnum_(-1) {}
        Cpoint(int f) : dist_(-1.0), density_(0), fnum_(f), nidx_(-1), oidx_(-1), cnum_(-1) {}
        Cpoint(Cpoint const& rhs) : dist_(rhs.dist_), density_(rhs.density_), fnum_(rhs.fnum_),
                                    nidx_(rhs.nidx_), oidx_(rhs.oidx_), cnum_(rhs.cnum_) {}
        Cpoint& operator=(Cpoint const& rhs) {
          if (&rhs != this) {
            dist_ = rhs.dist_; density_ = rhs.density_; fnum_ = rhs.fnum_;
            nidx_ = rhs.nidx_; oidx_ = rhs.oidx_; cnum_ = rhs.cnum_;
          }
          return *this;
        }
        /// Used to sort Carray by density
        //bool operator<(Cpoint const& second) const
        struct density_sort {
          inline bool operator()(Cpoint const& first, Cpoint const& second) const {
            if (first.density_ == second.density_)
              return (first.fnum_ < second.fnum_);
            else
              return (first.density_ < second.density_);
          }
        };
        /// Sort by density, descending
//        struct density_sort_descend {
//          inline bool operator()(Cpoint const& first, Cpoint const& second) const {
//              return (first.density_ > second.density_);
//          }
//        };
        /// Used to sort Carray by cluster number
        struct cnum_sort {
          inline bool operator()(Cpoint const& first, Cpoint const& second) const {
            if (first.cnum_ == second.cnum_)
              return (first.fnum_ < second.fnum_);
            else
              return (first.cnum_ < second.cnum_);
          }
        };
        /// Used to sort Carray by distance
        struct dist_sort {
          inline bool operator()(Cpoint const& first, Cpoint const& second) const {
            return (first.dist_ < second.dist_);
          }
        };
        double Dist()    const { return dist_; }
        int Density()    const { return density_; }
        //double Density() const { return density_; }
        int Fnum()       const { return fnum_; }
        int NearestIdx() const { return nidx_; }
        int Oidx()       const { return oidx_; }
        int Cnum()       const { return cnum_; }
        void SetDensity(int d)    { density_ = d; }
        void SetDist(double d)    { dist_ = d; }
        void SetNearestIdx(int n) { nidx_ = n; }
        void SetCluster(int c)    { cnum_ = c; }
        void AddDensity(double d) { density_ += d; }
      private:
        double dist_; ///< minimum distance to point with higher density
        int density_; ///< # other points within epsilon
        //double density_;
        int fnum_;    ///< Frame number.
        int nidx_;    ///< Index in Carray of nearest neighbor with higher density.
        int oidx_;    ///< Original index in Carray before sorting.
        int cnum_;    ///< Cluster number. -1 is no cluster.
    };
    typedef std::vector<Cpoint> Carray;
    Carray Points_; ///< Hold info for each point to be clustered.
};
#endif
