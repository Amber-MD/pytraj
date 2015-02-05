#ifndef INC_MULTIDIHEDRAL_H
#define INC_MULTIDIHEDRAL_H
#include "Action.h"
#include "DihedralSearch.h"
class Action_MultiDihedral : public Action {
  public:
    Action_MultiDihedral();
    static DispatchObject* Alloc() { return (DispatchObject*)new Action_MultiDihedral(); }
    static void Help();
  private:
    Action::RetType Init(ArgList&, TopologyList*, DataSetList*, DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);
    void Print();

    int debug_;
    DihedralSearch dihSearch_;    ///< Used to search for specified dihedrals
    std::vector<DataSet*> data_;  ///< Output DataSets, 1 per dihedral.
    bool range360_;
    Range resRange_;              ///< Residues to search for dihedrals.
    std::string dsetname_;
    DataFile* outfile_;
    // TODO: Replace these with new DataSet type
    DataSetList* masterDSL_;
};
#endif 
