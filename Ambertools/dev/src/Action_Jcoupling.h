#ifndef INC_ACTION_JCOUPLING_H
#define INC_ACTION_JCOUPLING_H
#include <vector>
#include <map>
#include <string>
#include "Action.h"
// Class: Action_Jcoupling
/// Calculate j-coupling values based on dihedrals using the Karplus equation.
/*! \author Original Code: J. Maier, Stony Brook 2011.
    \author Adapted by: D. Roe, Rutgers 2011.

    Two types of J-couplings are calculated by this code. If the type is
    0, the form used is that described by Chou et al. JACS (2003) 125 
    p.8959-8966 and the four constants represent A, B, C, and D. If the
    type is 1 (denoted C in $AMBERHOME/Karplus.txt) the form use is that 
    described by Perez et al. JACS (2001) 123 p.7081-7093 and the first 
    three constants represent C0, C1, and C2.
 */
class Action_Jcoupling: public Action {
  public:
    Action_Jcoupling();
    ~Action_Jcoupling();
    static DispatchObject* Alloc() { return (DispatchObject*)new Action_Jcoupling(); }
    static void Help();
  private:
    Action::RetType Init(ArgList&, TopologyList*, DataSetList*, DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);
    void Print() {}

    /// Load Karplus parameters from a file
    int loadKarplus(std::string);

    /// Hold Karplus parameters for a dihedral
    struct karplusConstant {
      NameType atomName[4]; ///< Name of each atom involved in dihedral
      int offset[4];        ///< Offsets
      double C[4];          ///< Constants
      int type;             ///< Calculation type (0=Chou, 1=Perez)
    };
    /// Hold all Karplus constants for a given residue
    typedef std::vector<karplusConstant> karplusConstantList;
    /// Map residue names to Karplus constants
    typedef std::map<std::string,karplusConstantList*> karplusConstantMap;
    /// Hold Karplus constants for all residues
    karplusConstantMap KarplusConstants_;
    /// Hold info for single j-coupling calculation
    struct jcouplingInfo {
      int residue; ///< Residue number
      int atom[4]; ///< Atom #s of the dihedral
      double *C;   ///< Pointer to C in associated karplusConstant structure
      int type;    ///< Calculation type (0=Chou, 1=Perez)
      DataSet* data_; ///< Hold Jcoupling vs frame
    };
    /// Hold info for all j-coupling calcs
    std::vector<jcouplingInfo> JcouplingInfo_;

    AtomMask Mask1_;         ///< Mask to search for dihedrals in.
    int debug_;              ///< Debug level.
    int Nconstants_;
    Topology* CurrentParm_;
    CpptrajFile outputfile_;
    DataFile* outfile_;
    // TODO: Replace these with new DataSet type
    DataSetList* masterDSL_;
    std::string setname_;
    int setcount_;
};
#endif  
