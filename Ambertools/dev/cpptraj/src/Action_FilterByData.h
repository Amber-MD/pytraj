#ifndef INC_ACTION_FILTERBYDATA_H
#define INC_ACTION_FILTERBYDATA_H
#include "Action.h"
#include "Array1D.h"
/// Filter out frames by DataSet
class Action_FilterByData : public Action {
  public:
    Action_FilterByData() : maxmin_(0) {}
    static DispatchObject* Alloc() { return (DispatchObject*)new Action_FilterByData(); }
    static void Help();
    /// For running as a separate command.
    size_t DetermineFrames() const;
    Action::RetType Init(ArgList&, TopologyList*, DataSetList*, DataFileList*, int);
    Action::RetType DoAction(int, Frame*, Frame**);
  private:
    Action::RetType Setup(Topology*, Topology**) { return Action::OK; }
    void Print() {}

    std::vector<double> Max_;
    std::vector<double> Min_;
    Array1D Dsets_;
    DataSet* maxmin_;
};
#endif
