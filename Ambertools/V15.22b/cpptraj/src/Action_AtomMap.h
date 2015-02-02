#ifndef INC_ACTION_ATOMMAP_H
#define INC_ACTION_ATOMMAP_H
#include "Action.h"
#include "AtomMap.h"
// Class: Action_AtomMap
/// Action used to map one molecule to another using AtomMaps
class Action_AtomMap : public Action {
  public:
    Action_AtomMap(); 
    static DispatchObject* Alloc() { return (DispatchObject*)new Action_AtomMap(); }
    static void Help();
    ~Action_AtomMap();
  private:
    Action::RetType Init(ArgList&, TopologyList*, FrameList*, DataSetList*,
                          DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);
    void Print() {}

    int debug_;
    AtomMap RefMap_;
    ReferenceFrame RefFrame_;
    AtomMap TgtMap_;
    ReferenceFrame TgtFrame_;

    std::vector<int> AMap_;
    bool maponly_;
    Frame* newFrame_;
    Topology* newParm_;
    Topology* stripParm_; // For stripping reference

    Frame rmsRefFrame_;
    Frame rmsTgtFrame_;
    bool rmsfit_;
    DataSet* rmsdata_;

    int mapBondsToUnique(AtomMap&, AtomMap&);
    int mapChiral(AtomMap&, AtomMap&);
    int mapByIndex(AtomMap&, AtomMap&);
    int mapUniqueRefToTgt(AtomMap&, AtomMap&, int);
    int MapAtoms(AtomMap&, AtomMap&);
    int MapUniqueAtoms(AtomMap&, AtomMap&);
    int MapWithNoUniqueAtoms( AtomMap&, AtomMap& );
};
#endif
