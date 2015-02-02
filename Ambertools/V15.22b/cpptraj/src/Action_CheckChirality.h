#ifndef INC_ACTION_CHECKCHIRALITY_H
#define INC_ACTION_CHECKCHIRALITY_H
#include "Action.h"
#include "Array1D.h"
/// Determine if amino acids are D or L 
class Action_CheckChirality: public Action {
  public:
    Action_CheckChirality();
    static DispatchObject* Alloc() { return (DispatchObject*)new Action_CheckChirality(); }
    static void Help();
  private:
    Action::RetType Init(ArgList&, TopologyList*, FrameList*, DataSetList*,
                          DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);
    void Print();

    struct ResidueInfo {
//      DataSet* data_; ///< data
      int num_; ///< residue number
      int isActive_;
      int n_;   ///< N coord index
      int ca_;  ///< C alpha coord index
      int c_;   ///< C coord index
      int cb_;  ///< C beta coord index
      int N_L_; ///< # times this residue was L
      int N_D_; ///< # times this residue was D
    };

    typedef std::vector<ResidueInfo> Rarray;
    Rarray resInfo_; 
    AtomMask Mask1_;
    std::string outfilename_;
    int ensembleNum_;
//    DataFile* outfile_;
//    std::string setname_;
//    DataSetList* masterDSL_;
};
#endif
