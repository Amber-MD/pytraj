#ifndef INC_ACTION_BOUNDS_H
#define INC_ACTION_BOUNDS_H
#include "Action.h"
/// Report the min/max XYZ values for atoms in mask.
class Action_Bounds : public Action {
  public:
    Action_Bounds();
    static DispatchObject* Alloc() { return (DispatchObject*)new Action_Bounds(); }
    static void Help();
  private:
    Action::RetType Init(ArgList&, TopologyList*, DataSetList*, DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);
    void Print();
    AtomMask mask_;
    std::string outfilename_;
    double max_[3];
    double min_[3];
    Vec3 dxyz_;
    int ensembleNum_;
    int offset_;
    DataSet* grid_;
};
#endif
