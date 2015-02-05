#ifndef INC_SYMMETRICRMSD_H
#define INC_SYMMETRICRMSD_H
#include "Action.h"
#include "ReferenceAction.h"
#include "SymmetricRmsdCalc.h"
/// Action to calculate symmetry-corrected RMSD
class Action_SymmetricRmsd : public Action {
  public:
    Action_SymmetricRmsd();
    static DispatchObject* Alloc() { return (DispatchObject*)new Action_SymmetricRmsd(); }
    static void Help();
  private:
    Action::RetType Init(ArgList&, TopologyList*, DataSetList*, DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);
    void Print() {}

    DataSet* rmsd_;           ///< Output DataSet
    ReferenceAction REF_;     ///< Hold reference frame/traj/options
    SymmetricRmsdCalc SRMSD_; ///< Symmetric RMSD calculation.
    AtomMask tgtMask_;        ///< Atom mask selecting target atoms.
    Frame selectedTgt_;       ///< Frame holding selected target atoms.
    bool remap_;              ///< If true, re-map target frame.
    Frame remapFrame_;        ///< Original target frame re-mapped for symmetry
    SymmetricRmsdCalc::Iarray targetMap_; ///< targetMap_[oldTgt] = newTgt
};
#endif
