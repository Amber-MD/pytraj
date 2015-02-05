#ifndef INC_ACTION_NASTRUCT_H
#define INC_ACTION_NASTRUCT_H
// Action_NAstruct
#include <vector>
#include <map>
#include "Action.h"
#include "AxisType.h"
#include "Range.h"
#include "DataSet_1D.h"
/// Basic Nucleic acid structure analysis. 
/** Calculate nucleic acid base/base pair structural parameters.
  * Algorithms for calculation of base/base pair structural parameters
  * adapted from:
  *   Babcock MS, Pednault EPD, Olson WK, "Nucleic Acid Structure Analysis: 
  *   Mathematics for Local Cartesian and Helical Structure Parameters That
  *   Are Truly Comparable Between Structures", J. Mol. Biol. (1994) 237,
  *   125-156.
  * NA base reference frame coordinates taken from:
  *   Olson WK, Bansal M, Burley SK, Dickerson RE, Gerstein M, Harvey SC,
  *   Heinemann U, Lu XJ, Neidle S, Shekked Z, Sklenar H, Suzuki M, Tung CS,
  *   Westhof E, Wolberger C, Berman H, "A Standard Reference Frame for the 
  *   Description of Nucleic Acid Base-pair Geometry", J. Mol. Biol. (2001)
  *   313, 229-237.
  */
class Action_NAstruct: public Action {
  public:
    Action_NAstruct();
    ~Action_NAstruct();
    static DispatchObject* Alloc() { return (DispatchObject*)new Action_NAstruct(); }
    static void Help();
  private:
    Action::RetType Init(ArgList&, TopologyList*, DataSetList*, DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);
    void Print();
    // Functions
    void ClearLists();

    int setupBaseAxes(Frame const&);

    int GCpair(NA_Base const&, NA_Base const&);
    int ATpair(NA_Base const&, NA_Base const&);
    int CalcNumHB(NA_Base const&, NA_Base const&);
    int determineBasePairing();

    int calculateParameters(NA_Axis const&, NA_Axis const&, NA_Axis*, double*);
    int helicalParameters(NA_Axis const&, NA_Axis const&, double *);
    int determineBaseParameters(int);
    void CalcPucker(NA_Base const&, int, int);
    int determineBasepairParameters(int);

    // Variables
    enum PmethodType { ALTONA=0, CREMER };
    PmethodType puckerMethod_;
    std::vector<NA_Base> Bases_;        ///< Hold ref/input coordinates for each base
    std::vector<NA_Axis> BaseAxes_;     ///< Hold axis coordinates for each base
    std::vector<NA_Axis> BasePairAxes_; ///< Hold axis coordinates for each base pair
    std::vector<int> NumberOfHbonds_;   ///< # hbonds between each base pair
    double HBdistCut2_;                 ///< distance Cutoff^2 for determining hydrogen bonds
    double originCut2_;                 ///< Cutoff^2 for determining base-pairing vi origins
    int maxResSize_;                    ///< Max residue size, used to set up frames for RMS fit.
    int debug_;
    int ensembleNum_;
    Range resRange_;                    ///< Range to search for NA residues.
    bool printheader_;                  ///< If true, print header to naout files.
    bool useReference_;                 ///< If true, use reference to determine base pairing.
    std::string outputsuffix_;          ///< Output file suffix (BP.<suffix> etc)
    std::string dataname_;              ///< NA DataSet name (default NA).

    typedef std::map<std::string, NA_Base::NAType> ResMapType;
    ResMapType CustomMap_;
    // Datasets - 1 entry per BP/BPstep
    typedef std::vector<DataSet_1D*> Darray;
    Darray PUCKER_;
    Darray SHEAR_;
    Darray STRETCH_;
    Darray STAGGER_;
    Darray BUCKLE_;
    Darray PROPELLER_;
    Darray OPENING_;
    Darray BPHBONDS_;
    Darray MAJOR_;
    Darray MINOR_;
    Darray SHIFT_;
    Darray SLIDE_;
    Darray RISE_;
    Darray TILT_;
    Darray ROLL_;
    Darray TWIST_;
    Darray XDISP_;
    Darray YDISP_;
    Darray HRISE_;
    Darray INCL_;
    Darray TIP_;
    Darray HTWIST_;
    // TODO: Replace these with new DataSet type
    DataSetList* masterDSL_;
#   ifdef NASTRUCTDEBUG
    // DEBUG - used to trigger AxisPDBwriter for first call of calculateParameters
    bool calcparam_;
#   endif
};
#endif
