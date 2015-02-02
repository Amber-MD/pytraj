#ifndef INC_ACTION_DISTRMSD_H
#define INC_ACTION_DISTRMSD_H
#include "Action.h"
#include "ReferenceAction.h"
// Class: Action_DistRmsd
/// Action to calculate the distance RMSD between frame and a reference frame.
class Action_DistRmsd: public Action, ReferenceAction {
  public:
    Action_DistRmsd();
    static DispatchObject* Alloc() { return (DispatchObject*)new Action_DistRmsd(); }
    static void Help();
  private:
    Action::RetType Init(ArgList&, TopologyList*, FrameList*, DataSetList*,
                          DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);
    void Print() {}

    DataSet *drmsd_;    ///< DRMSD DataSet
    AtomMask TgtMask_;  ///< Target mask.
    Frame SelectedTgt_; ///< Hold only target coords selected by TgtMask
};
#endif
