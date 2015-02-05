#ifndef INC_ACTION_RADGYR_H
#define INC_ACTION_RADGYR_H
#include "Action.h"
// Class: Action_Radgyr
/// Action to calculate the radius of gyration of atoms within a mask.
class Action_Radgyr: public Action {
  public:
    Action_Radgyr();
    static DispatchObject* Alloc() { return (DispatchObject*)new Action_Radgyr(); }
    static void Help();
  private:
    Action::RetType Init(ArgList&, TopologyList*, DataSetList*, DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);
    void Print() {}

    DataSet* rog_;
    DataSet* rogmax_;
    DataSet* rogtensor_;
    AtomMask Mask1_;
    bool calcRogmax_;
    bool calcTensor_;
    bool useMass_;
};
#endif
