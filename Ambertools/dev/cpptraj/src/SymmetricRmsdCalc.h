#ifndef INC_SYMMETRICRMSDCALC_H
#define INC_SYMMETRICRMSDCALC_H
#include "Topology.h"
#include "Hungarian.h"
#include "AtomMap.h"
/// Class for performing symmetry-corrected RMSD calculations.
class SymmetricRmsdCalc {
  public:
    typedef AtomMap::Iarray Iarray;
    SymmetricRmsdCalc();
    SymmetricRmsdCalc(AtomMask const&, bool, bool, Topology const&, int);
    /// Set fit, mass, and debug options.
    int InitSymmRMSD(bool, bool, int);
    /// Find symmetric atoms given a Topology and AtomMask.
    int SetupSymmRMSD(Topology const&, AtomMask const&, bool);
    /// Calc symm. RMSD using target and reference containing selected atoms only. 
    double SymmRMSD(Frame const&, Frame&);
    /// Calc symm. RMSD using target and pre-centered reference containing selected atoms.
    double SymmRMSD_CenteredRef(Frame const&, Frame const&);
    bool Fit()                    const { return fit_;         }
    bool UseMass()                const { return useMass_;     }
    Matrix_3x3 const& RotMatrix() const { return rotMatrix_;   }
    Vec3 const& TgtTrans()        const { return tgtTrans_;    }
    Iarray const& AMap()          const { return AMap_;        }
  private:
    typedef std::vector<Iarray> AtomIndexArray;
    /// Array of groups of potentially symmetric atoms
    AtomIndexArray SymmetricAtomIndices_;
    int debug_;
    Hungarian cost_matrix_; ///< Hungarian algorithm cost matrix.
    Iarray AMap_;           ///< AMap_[oldSelectedTgt] = newSelectedTgt
    Frame tgtRemap_;        ///< Selected target atoms re-mapped for symmetry.
    Matrix_3x3 rotMatrix_;  ///< Hold best-fit rotation matrix for target.
    Vec3 tgtTrans_;         ///< Hold translation of target to origin.
    bool fit_;              ///< If true, perform RMS best-fit.
    bool useMass_;          ///< If true, mass-weight calc.
};
#endif
