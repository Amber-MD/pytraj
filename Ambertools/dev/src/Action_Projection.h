#ifndef INC_ACTION_PROJECTION_H
#define INC_ACTION_PROJECTION_H
#include "Action.h"
#include "DataSet_Modes.h"
#include "ActionFrameCounter.h"
#include "Array1D.h"
/// project snapshots on normal modes
class Action_Projection : public Action, ActionFrameCounter {
  public:
    Action_Projection();
    static DispatchObject* Alloc() { return (DispatchObject*)new Action_Projection(); }
    static void Help();
  private:
    Action::RetType Init(ArgList&, TopologyList*, DataSetList*, DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);
    void Print() {}

    typedef std::vector<DataSet*> Darray;
    Darray project_;
    DataSet_Modes* modinfo_;
    int beg_;
    int end_;
    std::vector<double> sqrtmasses_;
    AtomMask mask_;
    Array1D DihedralSets_;
};
#endif
