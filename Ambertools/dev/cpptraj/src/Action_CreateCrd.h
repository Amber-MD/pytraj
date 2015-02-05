#ifndef ACTION_CREATECRD_H
#define ACTION_CREATECRD_H
#include "Action.h"
#include "DataSet_Coords.h"
class Action_CreateCrd : public Action {
  public:
    Action_CreateCrd();
    static DispatchObject* Alloc() { return (DispatchObject*)new Action_CreateCrd(); }
    static void Help();
  private:
    Action::RetType Init(ArgList&, TopologyList*, DataSetList*, DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);
    void Print() {}

    DataSet_Coords* coords_;
    int pindex_;
};
#endif
