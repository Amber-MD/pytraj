#ifndef INC_ACTION_CHANNEL_H
#define INC_ACTION_CHANNEL_H
#include "Action.h"
/// Experimental action for calculating solvent channels.
class Action_Channel : public Action {
  public:
    Action_Channel();
    static DispatchObject* Alloc() { return (DispatchObject*)new Action_Channel(); }
    static void Help();
  private:
    Action::RetType Init(ArgList&, TopologyList*, DataSetList*, DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);
    void Print() {}

    DataSet* grid_;
    AtomMask soluteMask_;
    AtomMask solventMask_;
    Vec3 dxyz_;
    std::vector<double> radii_;
};
#endif
