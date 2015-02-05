#ifndef INC_ACTION_STRIP_H
#define INC_ACTION_STRIP_H
#include "Action.h"
// Class: Action_Strip
/// Used to remove atoms from the state.
class Action_Strip: public Action {
  public:
    Action_Strip();
    static DispatchObject* Alloc() { return (DispatchObject*)new Action_Strip(); }
    static void Help();
    ~Action_Strip();
  private:
    Action::RetType Init(ArgList&, TopologyList*, DataSetList*, DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);
    void Print() {}

    Topology *oldParm_;
    Topology *newParm_;
    Frame newFrame_;
    std::string prefix_;
    AtomMask M1_;
    bool removeBoxInfo_;
};
// Class: Action_Unstrip
/// Signals to ActionList that the original traj parm should be restored.
class Action_Unstrip: public Action {
  public:
    Action_Unstrip() {}
    static DispatchObject* Alloc() { return (DispatchObject*)new Action_Unstrip(); }
    static void Help();

  private:
    Action::RetType Init(ArgList&, TopologyList*, DataSetList*, DataFileList*, int)
      { return Action::OK; }
    Action::RetType Setup(Topology*, Topology**)   { return Action::USEORIGINALFRAME; }
    Action::RetType DoAction(int, Frame*, Frame**) { return Action::USEORIGINALFRAME; }
    void Print() {}
};
#endif  
