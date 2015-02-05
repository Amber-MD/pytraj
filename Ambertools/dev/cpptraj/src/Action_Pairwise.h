#ifndef INC_ACTION_PAIRWISE_H
#define INC_ACTION_PAIRWISE_H
#include "Action.h"
#include "PDBfile.h"
#include "DataSet_MatrixDbl.h"
// Class: Pairwise 
/// Action to calculate nonbonded energy between pairs of atoms.
/** Functions in two ways:
  * - Calculate Lennard-Jones and Coulomb energy for each frame. 
  *   Also calculate the cumulative LJ and Coulomb energy on each atom.
  * - Calculate the Lennard-Jones and Coulomb energy between each
  *   pair of atoms in a reference structure. Calculate the difference
  *   in each pair from frame to reference (d = Ref - Frame). 
  */
class Action_Pairwise: public Action {
  public:
    Action_Pairwise();
    static DispatchObject* Alloc() { return (DispatchObject*)new Action_Pairwise(); }
    static void Help();
  private:
    Action::RetType Init(ArgList&, TopologyList*, DataSetList*, DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);
    void Print();

    typedef std::vector<double> Darray;
    enum EoutType {VDWOUT = 0, ELECOUT};
    enum PairCalcType { SET_REF, COMPARE_REF, NORMAL };
    PairCalcType nb_calcType_; ///< Type of nonbonded calc being performed
    AtomMask Mask0_;           ///< Calculate energy for atoms in mask
    AtomMask RefMask_;         ///< Reference mask
    Topology* CurrentParm_;    ///< Set to the current topology file.
    int N_ref_interactions_;   ///< Number of interactions in Ref w/exclusions
    int nframes_;              ///< Number of frames
    DataSet* ds_vdw_;          ///< Evdw dataset
    DataSet* ds_elec_;         ///< Eelec dataset
    DataSet_MatrixDbl* vdwMat_; ///< van der Waals energy map
    DataSet_MatrixDbl* eleMat_; ///< Coulomb energy map
    double ELJ_;               ///< Total VDW energy over all selected atoms.
    double Eelec_;             ///< Total elec. energy over all selected atoms.
    double cut_evdw_;          ///< van der Waals energy cutoff
    Darray atom_evdw_;         ///< Cumulative Evdw on each atom
    double cut_eelec_;         ///< Coulomb energy cutoff
    Darray atom_eelec_;        ///< Cumulative Eelec on each atom
    std::string mol2Prefix_;   ///< Mol2 file prefix for atoms satisfying cutoffs
    std::string avgout_;       ///< Filename for printing final results
    PDBfile PdbOut_;           ///< PDB with atoms colored by evdw/eelec
    CpptrajFile Eout_;         ///< Output file for atom energies.
    static const double QFAC;  ///< Convert charges to kcal/mol units
    /// Hold nonbond energy for a given atom pair
    struct NonbondEnergyType {
      double evdw;
      double eelec;
    };
    /// Hold nonbond energy for each pair of atoms in reference
    std::vector<NonbondEnergyType> ref_nonbondEnergy_;
    /// Hold cumulative LJ and elec energy for each atom
    //std::vector<NonbondEnergyType> atom_nonbondEnergy;

    /// Set up nonbondParm for given Parm and atoms in mask
    int SetupNonbondParm(AtomMask const&, Topology const&);
    /// Write energies to file
    inline void WriteEnergies(Topology const&, int, int, double, double, const char*);
    /// Calculate nonbond energy using nonbondParm for given frame
    void NonbondEnergy(Frame const&, Topology const&, AtomMask const&);
    /// Write mol2 file with atoms satisfying cutoff
    int WriteCutFrame(int, Topology const&, AtomMask const&, Darray const&, 
                      Frame const&, std::string const&);
    /// Write atoms satisfying cutoff for given energy type
    int PrintCutAtoms(Frame const&, int, EoutType, Darray const&, double);
};
#endif  
