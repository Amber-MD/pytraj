#ifndef INC_ACTION_AVERAGE_H
#define INC_ACTION_AVERAGE_H
#include "Action.h"
#include "ActionFrameCounter.h"
// Class: Action_Average
/// Sum up all coordinates and print the averaged coords in given format.
class Action_Average: public Action, ActionFrameCounter {
  public:
    Action_Average();
    static DispatchObject* Alloc() { return (DispatchObject*)new Action_Average(); }
    static void Help();
    ~Action_Average();
  private:
    Action::RetType Init(ArgList&, TopologyList*, DataSetList*, DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);
    void Print();

    int ensembleNum_;
    int debug_;
    AtomMask Mask1_;
    Frame* AvgFrame_;
    Topology AvgParm_;
    ArgList trajArgs_;
    int Natom_;
    int Nframes_;
    std::string avgfilename_;
    DataSet* crdset_;         ///< DataSet to save avg coords to.
};
#endif  
