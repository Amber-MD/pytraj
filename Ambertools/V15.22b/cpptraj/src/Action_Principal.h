#ifndef INC_ACTION_PRINCIPAL_H
#define INC_ACTION_PRINCIPAL_H
#include "Action.h"
class Action_Principal : public Action {
  public:
    Action_Principal();
    static DispatchObject* Alloc() { return (DispatchObject*)new Action_Principal(); }
    static void Help();
  private:
    bool doRotation_;
    bool useMass_;
    int debug_;
    AtomMask mask_;
    CpptrajFile outfile_;

    Action::RetType Init(ArgList&, TopologyList*, FrameList*, DataSetList*,
                          DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);
    void Print() {}
};
#endif
