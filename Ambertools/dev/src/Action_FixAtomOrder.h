#ifndef INC_ACTION_FIXATOMORDER_H
#define INC_ACTION_FIXATOMORDER_H
#include "Action.h"
// Class: Action_FixAtomOrder
/// Fix atom ordering in parm where atoms in mols are not sequential. 
class Action_FixAtomOrder: public Action {
  public:
    Action_FixAtomOrder();
    ~Action_FixAtomOrder();
    static DispatchObject* Alloc() { return (DispatchObject*)new Action_FixAtomOrder(); }
    static void Help();
  private:
    Action::RetType Init(ArgList&, TopologyList*, DataSetList*, DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);
    void Print() {}

    void VisitAtom(int,int,Topology const&);

    int debug_;
    typedef std::vector<int> MapType;
    MapType atomMap_;
    MapType molNums_; ///< Hold molecule number for each atom.
    Topology* newParm_;
    Frame newFrame_;
    std::string prefix_;
};
#endif
