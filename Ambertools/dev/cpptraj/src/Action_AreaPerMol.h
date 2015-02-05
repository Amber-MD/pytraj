#ifndef INC_ACTION_AREAPERMOL_H
#define INC_ACTION_AREAPERMOL_H
#include "Action.h"
/// Calculate area per molecule for given molecule/dimensions. 
class Action_AreaPerMol: public Action {
  public:
    Action_AreaPerMol();
    static DispatchObject* Alloc() { return (DispatchObject*)new Action_AreaPerMol(); }
    static void Help();
  private:
    enum AreaType { XY, XZ, YZ };

    Action::RetType Init(ArgList&, TopologyList*, DataSetList*, DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);
    void Print() {}

    DataSet *area_per_mol_;
    double Nmols_;
    double Nlayers_;
    AreaType areaType_;
    AtomMask Mask1_;
};
#endif
