#ifndef INC_ACTION_CLUSTERDIHEDRAL_H
#define INC_ACTION_CLUSTERDIHEDRAL_H
#include "Action.h"
class Action_ClusterDihedral : public Action {
  public:
    Action_ClusterDihedral();
    static DispatchObject* Alloc() { return (DispatchObject*)new Action_ClusterDihedral(); }
    static void Help();
  private:
    Action::RetType Init(ArgList&, TopologyList*, DataSetList*, DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);
    void Print();
    int ReadDihedrals(std::string const&);
    class DCnode;
    class DCmask;

    std::vector<DCnode> dcarray_;  ///< Hold counts for each bin# combo.
    std::vector<DCmask> DCmasks_;  ///< Hold 4 atom mask for each dihedral
    std::vector<int> Bins_;        ///< For holding bin combo each frame
    int phibins_;
    int psibins_;
    int CUT_;
    int lastframe_;
    Topology* dcparm_;             ///< Set to first parm encountered.
    std::string outfile_;
    std::string framefile_; // filenames[1]
    std::string infofile_;  // filenames[2]
    AtomMask mask_;
    DataSet* CVT_; ///< Hold # clusters vs time.
    double minimum_; ///< Value of first bin 
    int debug_;
    int ensembleNum_;
};
// -----------------------------------------------------------------------------
class Action_ClusterDihedral::DCnode {
  public:
    DCnode() : count_(0) {}
    DCnode(std::vector<int>& binIn, int frameIn) : 
           BinIDs_(binIn), frames_(1, frameIn), count_(1) {}
    DCnode(const DCnode& rhs) : 
           BinIDs_(rhs.BinIDs_), frames_(rhs.frames_), count_(rhs.count_) {}
    DCnode& operator=(const DCnode& rhs) {
      if (this==&rhs) return *this;
      BinIDs_ = rhs.BinIDs_;
      frames_ = rhs.frames_;
      count_ = rhs.count_;
      return *this;
    }
    void Increment()        { ++count_;                 }
    void AddFrame(int fIn)  { frames_.push_back( fIn ); }
    // Want sort in descending order, so reverse '>'
    bool operator<(const DCnode& rhs) const  { return (count_ > rhs.count_);  }
    bool operator>(const DCnode& rhs) const  { return (count_ < rhs.count_);  }
    bool operator==(const DCnode& rhs) const { return (count_ == rhs.count_); }
    bool BinMatch(std::vector<int>& binIn) {
      return (std::equal(BinIDs_.begin(), BinIDs_.end(), binIn.begin()));
    }
    long int Count()      { return count_; }
    typedef std::vector<int>::const_iterator bin_it;
    bin_it binbegin()     { return BinIDs_.begin(); }
    bin_it binend()       { return BinIDs_.end();   }
    typedef std::vector<int>::const_iterator frame_it;
    frame_it framebegin() { return frames_.begin(); }
    frame_it frameend()   { return frames_.end();   }
    int NumFrames()       { return (int)frames_.size(); }
  private:
    std::vector<int> BinIDs_;
    std::vector<int> frames_;
    long int count_;
};
// -----------------------------------------------------------------------------
class Action_ClusterDihedral::DCmask {
  public: 
    DCmask() : a1_(0), a2_(0), a3_(0), a4_(0), bins_(0) {}
    DCmask(int a1, int a2, int a3, int a4, int bins, double min) :
           a1_(a1), a2_(a2), a3_(a3), a4_(a4), 
           bins_(bins), step_(360/(double)bins), min_(min) {}
    int A1()   { return a1_; }
    int A2()   { return a2_; }
    int A3()   { return a3_; }
    int A4()   { return a4_; }
    int Bins() { return bins_; }
    double Step() { return step_; }
    double Min() { return min_; }
  private:
    int a1_;
    int a2_;
    int a3_;
    int a4_;
    int bins_;
    double step_; 
    double min_;
};
#endif
