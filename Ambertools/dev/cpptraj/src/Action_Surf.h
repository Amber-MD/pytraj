#ifndef INC_ACTION_SURF_H
#define INC_ACTION_SURF_H
#include "Action.h"
// Class: Action_Surf
/// Calculate LCPO surface area.
/** LCPO method from:
  * -  J. Weiser, P.S. Shenkin, and W.C. Still,
  *    "Approximate atomic surfaces from linear combinations of pairwise
  *    overlaps (LCPO)", J. Comp. Chem. 20:217 (1999).
  */
class Action_Surf: public Action {
  public:
    Action_Surf();
    static DispatchObject* Alloc() { return (DispatchObject*)new Action_Surf(); }
    static void Help();
  private:
    Action::RetType Init(ArgList&, TopologyList*, DataSetList*, DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);
    void Print() {}

    DataSet* surf_;
    AtomMask Mask1_;
    AtomMask atomi_neighborMask_;
    AtomMask atomi_noNeighborMask_;
    AtomMask atomj_neighborMask_;
    /// Contain data for an atoms LCPO SA calc
    // TODO: Rework VDW storage
    struct SurfInfo {
      double vdwradii;
      double P1;
      double P2;
      double P3;
      double P4;
    };
    /// Contain LCPO data for all atoms in atomi_neighborMask
    std::vector<SurfInfo> SurfaceInfo_neighbor_;
    /// Contain LCPO data for all atoms in atomi_noNeighborMask
    std::vector<SurfInfo> SurfaceInfo_noNeighbor_;
    /// Contain vdw radii for all atoms
    std::vector<double> VDW_;

    void AssignLCPO(SurfInfo *, double, double, double, double , double );
    void SetAtomLCPO(Topology const&,int, SurfInfo*);
};
#endif
