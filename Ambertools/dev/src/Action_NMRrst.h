#ifndef INC_ACTION_NMRRST_H
#define INC_ACTION_NMRRST_H
#include "Action.h"
#include "ImagedAction.h"
#include "BufferedLine.h"
// Class: Action_NMRrst
class Action_NMRrst: public Action {
  public:
    Action_NMRrst();
    static DispatchObject* Alloc() { return (DispatchObject*)new Action_NMRrst(); }
    static void Help();
  private:
    Action::RetType Init(ArgList&, TopologyList*, DataSetList*, DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);
    void Print();

    int ReadNmrRestraints( std::string const& );
    int ReadXplor( BufferedLine& );
    int ReadAmber( BufferedLine& );
    /// Type to hold NOE data.
    struct noeDataType {
      int resNum1_;     ///< First residue number.
      int resNum2_;     ///< Second residue number.
      std::string aName1_; ///< First atom name.
      std::string aName2_; ///< Second atom name.
      AtomMask dMask1_; ///< First mask for distance pair.
      AtomMask dMask2_; ///< Second mask for distance pair.
      double bound_;    ///< Lower bound.
      double boundh_;   ///< Upper bound.
      double rexp_;     ///< Expected distance
      DataSet* dist_;   ///< Distance DataSet.
      bool active_;     ///< True if NOE was properly set up.
    };
    typedef std::vector<noeDataType> noeDataArray;
    noeDataArray NOEs_;

    typedef std::vector<int> Iarray;

    class Site;
    typedef std::vector<Site> SiteArray;

    /// Used to map NOEs to unique values, res1 always < res2. 
    typedef std::pair<int,int> Ptype;

    class NOEtype;
    typedef std::vector<NOEtype> NOEtypeArray;
    NOEtypeArray noeArray_;
    
    ImagedAction Image_;
    std::string setname_;
    std::string findOutputName_;
    AtomMask Mask_;
    DataSetList* masterDSL_; // TODO: Replace these with new DataSet type
    size_t numNoePairs_; ///< Used to check if # of pairs has changed
    double max_cut_; ///< Min distance cutoff for NOE to be considered
    double strong_cut_;
    double medium_cut_;
    double weak_cut_;
    int resOffset_;
    int debug_;
    int ensembleNum_;
    int nframes_; ///< Total # of frames.
    bool useMass_;
    bool findNOEs_;
    bool series_; ///< If true save NOE distance data.
};
// ----- Associated Classes ----------------------------------------------------
/// Potential NOE site.
class Action_NMRrst::Site {
  public:
    Site() : resNum_(-1) {}
    Site(int r, Iarray const& i) : resNum_(r), indices_(i), shortestCount_(i.size(), 0) {}
    int ResNum()                   const { return resNum_;           }
    unsigned int Nindices()        const { return indices_.size();   }
    int Idx(unsigned int i)        const { return indices_[i];       }
    int Count(unsigned int i)      const { return shortestCount_[i]; }
    void Increment(int c)                { ++shortestCount_[c];      }
    std::string SiteLegend(Topology const&) const;
  private:
    int resNum_; ///< Site residue number.
    Iarray indices_; ///< Site atom indices.
    Iarray shortestCount_; ///< # times atom was part of shortest distance.
};
/// NOE between two sites.
class Action_NMRrst::NOEtype {
  public:
    NOEtype() : dist2_(0), r6_avg_(0.0) {}
    NOEtype(Site const& s1, Site const& s2, DataSet* d, std::string const& l) :
      site1_(s1), site2_(s2), legend_(l), dist2_(d), r6_avg_(0.0) {}
    Site const& Site1()    const { return site1_;    }
    Site const& Site2()    const { return site2_;    }
    double R6_Avg()        const { return r6_avg_;   }
    const char* legend()   const { return legend_.c_str(); }
    DataSet* Data()              { return dist2_;    }
    void SetR6Avg(double r6)     { r6_avg_ = r6;     }
    std::string PrintNOE() const;
    void UpdateNOE(int i, double d2, unsigned int c1, unsigned int c2) {
      if (dist2_ != 0) {
        float fval = (float)d2;
        dist2_->Add(i, &fval);
      }
      site1_.Increment( c1 );
      site2_.Increment( c2 );
      r6_avg_ += ( 1.0 / (d2 * d2 * d2) );
    }
    inline bool operator<(const NOEtype& rhs) const { return (r6_avg_ < rhs.r6_avg_); } 
  private:
    Site site1_;     ///< First site, lower resNum
    Site site2_;     ///< Second site, higher resNum
    std::string legend_;
    DataSet* dist2_; ///< Distance^2 data.
    double r6_avg_;  ///< Sum of r^-6 over all frames.
};
#endif
