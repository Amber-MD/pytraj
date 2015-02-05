#ifndef INC_MULTIVECTOR_H
#define INC_MULTIVECTOR_H
#include "Action.h"
#include "DataSet_Vector.h"
#include "Range.h"
class Action_MultiVector : public Action {
  public:
    Action_MultiVector();
    static DispatchObject* Alloc() { return (DispatchObject*)new Action_MultiVector(); }
    static void Help();
  private:
    Action::RetType Init(ArgList&, TopologyList*, DataSetList*, DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);
    void Print() {}

    int debug_;
    std::vector<DataSet_Vector*> data_;  ///< Output DataSets, 1 per vector.
    Range resRange_;                     ///< Residues to search for vectors.
    std::string dsetname_;
    NameType name1_;
    NameType name2_;
    std::vector<int> CrdIdx1_;
    std::vector<int> CrdIdx2_;
    DataFile* outfile_;
    // TODO: Replace these with new DataSet type
    DataSetList* masterDSL_;
    bool ired_;
};
#endif 
