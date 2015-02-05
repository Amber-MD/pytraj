#ifndef INC_ACTION_DSSP_H
#define INC_ACTION_DSSP_H
#include "Action.h"
/// Calculate protein secondary structure using DSSP algorithm.
class Action_DSSP : public Action {
  public:
    Action_DSSP();
    static DispatchObject* Alloc() { return (DispatchObject*)new Action_DSSP(); }
    static void Help();
  private:
    Action::RetType Init(ArgList&, TopologyList*, DataSetList*, DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);
    void Print();
    // Enum and static vars
    //enum SStype { ALPHA=0, ANTI, PARA, H3_10, HPI, TURN, BEND, NONE };
    enum SStype { NONE=0, PARA, ANTI, H3_10, ALPHA, HPI, TURN, BEND };
    static const int NSSTYPE;      ///< # of secondary structure types.
    static const double DSSP_fac;  ///< DSSP factor for calc. Hbond energy.
    static const char dssp_char[]; ///< DSSP 1 character SS names
    static const char* SSchar[];   ///< PTRAJ 1 character SS names
    static const char* SSname[];   ///< Full SS names
    /// Hold SS-related data for each residue
    struct SSres {
      std::vector<int> CO_HN_Hbond; ///< 1 if this res CO bonded to res X NH.
      DataSet* resDataSet;      ///< DataSet for SS assignment each frame.
      int SSprob[8];            ///< Hold count for each SS type
      SStype sstype;            ///< Assigned secondary structure
      int C;                    ///< Coord idx of BB carbon
      int O;                    ///< Coord idx of BB oxygen
      int N;                    ///< Coord idx of BB nitrogen
      int H;                    ///< Coord idx of BB hydrogen
      int CA;                   ///< Coord idx of BB alpha carbon
      bool isSelected;          ///< True if calculating SS for this residue.
      bool hasCO;               ///< True if both C and O atoms selected.
      bool hasNH;               ///< True if both N and H atoms selected. 
    };
    std::vector<SSres> SecStruct_; ///< Hold SS-related data for all residues
    // Class variables
    int ensembleNum_;
    int debug_;
    DataFile* outfile_;       ///< Output Data file
    DataFile* dsspFile_;      ///< Sum output file
    std::string dsetname_;    ///< DSSP data set name
    std::string assignout_;   ///< Assignment output file.
    AtomMask Mask_;           ///< Mask used to determine selected residues
    int Nres_;                ///< Current total # of residues
    unsigned int Nselected_;  ///< Current # residues selected.
    int Nframe_;              ///< # of frames, for calculating SS avg.
    bool printString_;        ///< If true print 1 char per residue indicating ss type
    // TODO: Replace these with new type of DataSet
    DataSetList* masterDSL_;
    DataSet* totalDS_[8];
    NameType BB_N_;
    NameType BB_H_;
    NameType BB_C_;
    NameType BB_O_;
    NameType BB_CA_;
    // Private fns
    inline int isBonded(int, int);
    inline void SSassign(int, int, SStype, bool);
    static inline bool HasPriority(SStype, SStype);
};
#endif
