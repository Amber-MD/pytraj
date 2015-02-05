#ifndef INC_ACTION_CLOSEST_H
#define INC_ACTION_CLOSEST_H
#include "Action.h"
#include "ImagedAction.h"
// Class: Action_Closest
/// Modify the state so that only the closest solvent molecules are kept.
class Action_Closest: public Action, ImagedAction {
  public:
    Action_Closest();
    static DispatchObject* Alloc() { return (DispatchObject*)new Action_Closest(); }
    static void Help();
    ~Action_Closest();
  private:
    Action::RetType Init(ArgList&, TopologyList*, DataSetList*, DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);
    void Print() {}

    DataFile *outFile_;     ///< Output file for data on closest molecules
    DataSet *framedata_;    ///< Frame number for each closest molecule.
    DataSet *moldata_;      ///< Mol# for each closest molecule.
    DataSet *distdata_;     ///< Closest distance of each molecule.
    DataSet *atomdata_;     ///< First atom of each closest molecule.
    int Nclosest_;          ///< Index into Closest molecule DataSets.
    std::string prefix_;    ///< Output topology prefix.
    int closestWaters_;     ///< Closest # of molecules to keep.
    bool firstAtom_;        ///< If true just calc based on molecule first atom.
    AtomMask stripMask_;    ///< Mask including all solute and closest molecules.
    AtomMask distanceMask_; ///< Mask of atoms to calculate distance from solvent to.
    Topology *newParm_;     ///< New topology with solute and closest molecules.
    int NsolventMolecules_; ///< # of solvent molecules in SolventMols.
    int debug_;
    Frame newFrame_;        ///< New frame with solute and kept waters.
    typedef std::vector<int> Iarray;
    /// Hold atom #s of kept solvent in new frame for placing into stripMask.
    Iarray keptWaterAtomNum_;
    /** The moldist structure is used in order to preserve the original
      * solvent molecule numbers after sorting. */
    struct MolDist {
      int mol;        ///< Original solvent molecule number (starts from 1).
      double D;       ///< Closest distance of solvent molecule to atoms in distanceMask.
      AtomMask mask;  ///< Original topology solvent molecule atom mask.
      Iarray solventAtoms; ///< Actual solvent atom #s to loop over.
    };
    /// Return true if the first molecule is closer than the second
    struct moldist_cmp {
      inline bool operator()(MolDist const& first, MolDist const& second) const {
        return (first.D < second.D);
      }
    };
    std::vector<MolDist> SolventMols_;
};
#endif  
