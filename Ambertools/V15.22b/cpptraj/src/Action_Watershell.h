#ifndef INC_ACTION_WATERSHELL_H
#define INC_ACTION_WATERSHELL_H
#include "Action.h"
#include "ImagedAction.h"
/// Calculate number of solvent residues in 1st/2nd solvation shell.
class Action_Watershell : public Action, ImagedAction {
  public:
    Action_Watershell();
#   ifdef _OPENMP
    ~Action_Watershell();
#   endif
    static DispatchObject* Alloc() { return (DispatchObject*)new Action_Watershell(); }
    static void Help();
  private:
    Action::RetType Init(ArgList&, TopologyList*, FrameList*, DataSetList*,
                          DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);
    void Print() { return; }

    AtomMask soluteMask_;
    AtomMask solventMask_;
    std::string solventmaskexpr_;
    double lowerCutoff_;
    double upperCutoff_;
    Topology* CurrentParm_;
    DataSet* lower_;
    DataSet* upper_;
    int numthreads_;
#   ifdef _OPENMP
    int** activeResidues_thread_;
    int NactiveResidues_;
#   else
    std::vector<int> activeResidues_;
#   endif
};
#endif
