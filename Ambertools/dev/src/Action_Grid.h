#ifndef INC_ACTION_GRID_H
#define INC_ACTION_GRID_H
#include "Action.h"
#include "DataSet_GridFlt.h"
#include "GridAction.h"
class Action_Grid : public Action, private GridAction {
  public:
    Action_Grid();
    static DispatchObject* Alloc() { return (DispatchObject*)new Action_Grid(); }
    static void Help();
  private:
    Action::RetType Init(ArgList&, TopologyList*, DataSetList*, DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);
    void Print();

    void PrintPDB(double);

    enum NormType { NONE=0, TO_FRAME, TO_DENSITY };
    NormType normalize_;
    int ensembleNum_;
    double density_;
    double max_;
    double madura_;
    double smooth_;
    unsigned int nframes_;
    bool invert_;
    AtomMask mask_;
    std::string pdbname_;
    DataSet_GridFlt* grid_;
};
#endif
