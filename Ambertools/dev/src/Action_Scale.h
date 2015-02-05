#ifndef INC_ACTION_SCALE_H
#define INC_ACTION_SCALE_H
#include "Action.h"
class Action_Scale : public Action {
  public:
    Action_Scale();
    static DispatchObject* Alloc() { return (DispatchObject*)new Action_Scale(); }
    static void Help();
  private:
    Action::RetType Init(ArgList&, TopologyList*, DataSetList*, DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);
    void Print() {}

    AtomMask mask_;
    double sx_;
    double sy_;
    double sz_;
};
#endif
