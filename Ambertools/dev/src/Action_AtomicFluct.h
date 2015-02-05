#ifndef INC_ACTION_ATOMICFLUCT_H
#define INC_ACTION_ATOMICFLUCT_H
#include "Action.h"
#include "ActionFrameCounter.h"
class Action_AtomicFluct : public Action, ActionFrameCounter {
  public :
    Action_AtomicFluct();
    static DispatchObject* Alloc() { return (DispatchObject*)new Action_AtomicFluct(); }
    static void Help();
  private :
    Action::RetType Init(ArgList&, TopologyList*, DataSetList*, DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);
    void Print();

    enum outputType { BYATOM = 0, BYRES, BYMASK };

    int ensembleNum_;
    Frame SumCoords_;         ///< Hold the average coordinates.
    Frame SumCoords2_;        ///< Hold the variance of coordinates.
    Frame Cross_;             ///< Hold cross-terms for calculating covariance.
    AtomMask Mask_;
    int sets_;
    bool bfactor_;
    bool calc_adp_;
    std::string adpoutname_;
    std::string outfilename_;
    std::string setname_;
    Topology *fluctParm_;
    outputType outtype_;
    DataSet* dataout_;
    DataFile* outfile_;
};
#endif
