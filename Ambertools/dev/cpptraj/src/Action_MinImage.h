#ifndef INC_ACTION_MINIMAGE_H
#define INC_ACTION_MINIMAGE_H
#include "Action.h"
#include "ImagedAction.h"
//#incl ude "PDBfile.h" // DEBUG
/// Action to calculate minimum non-self distance between atoms in two masks.
class Action_MinImage: public Action {
  public:
    Action_MinImage();
    static DispatchObject* Alloc() { return (DispatchObject*)new Action_MinImage(); }
    static void Help();
  private:
    Action::RetType Init(ArgList&, TopologyList*, DataSetList*, DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);
    void Print() {}

    double MinNonSelfDist2(Vec3 const&, Vec3 const&);

    ImagedAction image_;
    Matrix_3x3 ucell_, recip_;
    DataSet* dist_;      ///< Will hold DataSet of calculated distances.
    DataSet* atom1_;
    DataSet* atom2_;
    bool useMass_;       ///< If true, mass-weight distances.
    bool calcUsingMask_; ///< If true use center of masks
    AtomMask Mask1_;
    AtomMask Mask2_;
    std::vector<double> minDist_;
    std::vector<int> minAtom1_;
    std::vector<int> minAtom2_;
    //PDBfile pdbout_; // DEBUG
};
#endif
