#ifndef INC_ACTION_SETVELOCITY_H
#define INC_ACTION_SETVELOCITY_H
#include "Action.h"
#include "Random.h"
/// Calculate the temperature of parts of a system.
class Action_SetVelocity : public Action {
  public:
    Action_SetVelocity();
    static DispatchObject* Alloc() { return (DispatchObject*)new Action_SetVelocity(); }
    static void Help();
  private:
    Action::RetType Init(ArgList&, TopologyList*, DataSetList*, DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);
    void Print() {}

    AtomMask Mask_;
    std::vector<double> SD_;
    double tempi_;
    Random_Number RN_;
};
#endif
