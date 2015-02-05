#ifndef INC_ACTION_DIHEDRAL_H
#define INC_ACTION_DIHEDRAL_H
#include "Action.h"
/// Calculate dihedral in a Frame
class Action_Dihedral: public Action {
  public:
    Action_Dihedral();
    static DispatchObject* Alloc() { return (DispatchObject*)new Action_Dihedral(); }
    static void Help();
  private:
    Action::RetType Init(ArgList&, TopologyList*, DataSetList*, DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);
    void Print();

    DataSet* dih_;
    bool useMass_;
    bool range360_;
    AtomMask M1_;
    AtomMask M2_;
    AtomMask M3_;
    AtomMask M4_;
};
#endif  
