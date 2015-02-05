#ifndef INC_ACTION_CHECKSTRUCTURE_H
#define INC_ACTION_CHECKSTRUCTURE_H
#include "Action.h"
#include "ImagedAction.h"
// Class: Action_CheckStructure 
/// Action to check bond lengths and bad overlaps between non-bonded atoms 
class Action_CheckStructure: public Action, ImagedAction {
  public:
    Action_CheckStructure();
    static DispatchObject* Alloc() { return (DispatchObject*)new Action_CheckStructure(); }
    static void Help();
    ~Action_CheckStructure();
    // These are made public for use in other actions (e.g. Action_DihedralScan)
    Action::RetType Init(ArgList&, TopologyList*, DataSetList*, DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    int CheckFrame(int, Frame const&);
    void Print() {}
  private:
    Action::RetType DoAction(int, Frame*, Frame**);

    void SetupBondlist(BondArray const&, BondParmArray const&, AtomMask const&);
    /// Used to cache bond parameters
    struct bond_list {
      double req; ///< Hold (req + bondoffset)^2
      int atom1;
      int atom2;
    };
    typedef std::vector<bond_list> BondListType;
    BondListType bondL_;
    /// Sort first by atom1, then by atom2
    struct bond_list_cmp {
      inline bool operator()(bond_list const& first, bond_list const& second) const {
        if (first.atom1 < second.atom1) {
          return true;
        } else if (first.atom1 == second.atom1) {
          if (first.atom2 < second.atom2) return true;
        } 
        return false;
      }
    };
#   ifdef _OPENMP
    /// Hold position in bondL for each atom.
    std::vector<BondListType::const_iterator> BondsToAtomBegin_;
    enum ProblemType { NONE = 0, BOND, DISTANCE, BOTH };
    class Problem;
    std::vector<Problem>* problemIndices_;
    int numthreads_;
#   endif
    AtomMask Mask1_;
    double bondoffset_;
    double nonbondcut2_;
    bool bondcheck_;
    bool silent_;
    bool skipBadFrames_;
    CpptrajFile outfile_;
    Topology* CurrentParm_;
    int debug_;
};
#ifdef _OPENMP
class Action_CheckStructure::Problem {
  public:
  Problem() : Dist_(0.0), frameNum_(-1), type_(BOND), atom1_(-1), atom2_(-1) {}
  Problem(int f) : Dist_(0.0), frameNum_(f), type_(BOND), atom1_(-1), atom2_(-1) {}
  Problem(Problem const& rhs) : 
    Dist_(rhs.Dist_), frameNum_(rhs.frameNum_), type_(rhs.type_), 
    atom1_(rhs.atom1_), atom2_(rhs.atom2_) {}
  double Dist_;
  int frameNum_; 
  int type_; // 0: NONE, 1: BOND, 2: DISTANCE, 3: BOTH
  int atom1_;
  int atom2_;
};
#endif
#endif  
