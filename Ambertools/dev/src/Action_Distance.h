#ifndef INC_ACTION_DISTANCE_H
#define INC_ACTION_DISTANCE_H
#include "Action.h"
#include "ImagedAction.h"
// Class: Action_Distance
/// Action to calculate a distance between atoms in two masks.
class Action_Distance: public Action, ImagedAction {
  public:
    Action_Distance();
    static DispatchObject* Alloc() { return (DispatchObject*)new Action_Distance(); }
    static const char* NOE_Help;
    static void Help();
    static int NOE_Args(ArgList& argIn, double&, double&, double&);
  private:
    Action::RetType Init(ArgList&, TopologyList*, DataSetList*, DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);
    void Print() {}

    DataSet* dist_;  ///< Will hold DataSet of calculated distances.
    bool useMass_;   ///< If true, mass-weight distances.
    AtomMask Mask1_;
    AtomMask Mask2_;
};
#endif  
