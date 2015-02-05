#ifndef INC_ACTION_OUTTRAJ_H
#define INC_ACTION_OUTTRAJ_H
// Action_Outtraj
#include "Action.h"
#include "Trajout.h"
#include "DataSet_1D.h"
/// Write out a trajectory inside the ActionList
class Action_Outtraj: public Action {
  public:
    Action_Outtraj();
    static DispatchObject* Alloc() { return (DispatchObject*)new Action_Outtraj(); }
    static void Help();
  private:
    Action::RetType Init(ArgList&, TopologyList*, DataSetList*, DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);
    void Print();

    Trajout outtraj_;
    Topology* CurrentParm_;
    std::vector<double> Max_;
    std::vector<double> Min_;
    std::vector<DataSet_1D*> Dsets_;
};
#endif
