#ifndef INC_ACTION_SPAM_H
#define INC_ACTION_SPAM_H
#include "Action.h"
#include "ImagedAction.h"
#include "Vec3.h"
/*
SPAM is a water profiling technique developed by Guanglei Cui at
GlaxoSmithKline (GSK). The original implementation involved a set of specialized
Python scripts interacting with VMD (via the VolMap tool), numpy, NAMD (for the
SPAM energy calculations) and R (for the free energy calculation using a
specialized kernel density estimate). While that implementation demonstrated
proof of principle, simply re-ordering the trajectory for use with NAMD proved
to be a performance bottleneck because it was written in Python. SPAM was
rewritten from the ground up in cpptraj, significantly improving efficiency and
providing a simpler interface.

The original C++ implementation of SPAM in cpptraj was done by Jason Swails
while interning at GSK. This code was built as a patch on top of cpptraj v.12
and was rewritten by Jason Swails for the current cpptraj version.

 (C) 2012 - 2013
*/

// Class: Action_Spam
/** Action to calculate the Linear Interaction Energy (effectively the nonbonded 
  * energies between two different masks
  */
class Action_Spam: public Action, ImagedAction {
  public:
    Action_Spam();
    static DispatchObject* Alloc() { return (DispatchObject*) new Action_Spam(); }
    static void Help();
  private:
    int ensembleNum_;
    /** \brief Name of the solvent residues */
    std::string solvname_;
    /** \brief SPAM free energy of the bulk solvent */
    double bulk_;
    /** \brief Determines if we are running a pure water simulation
      * to derive bulk properties
      */
    bool purewater_;
    /** \brief Flag indicating whether or not solvent should be reordered */
    bool reorder_;
    /** \brief Non-bonded cutoff in Angstroms (squared) */
    double cut2_;
    /** \brief 1 / cut2_ */
    double onecut2_;
    /** \brief twice the cutoff (to test if boxes are big enough) */
    double doublecut_;
    /** \brief Name of the SPAM info file */
    std::string infoname_;
    /** \brief Mask for selecting individual solvent residues */
    AtomMask mask_;
    /** \brief File containing the summary of all SPAM statistics */
    std::string summaryfile_;
    /** \brief File containing all SPAM energies for each site */
    std::string datafile_;
    /** \brief Size of the water site. This is a full edge length or diameter */
    double site_size_;
    /** \brief A list of all omitted frames for each peak */
    std::vector< std::vector<int> > peakFrameData_;
    /** \brief The topology instance so we can extract necessary parameters for
      * energy evaluations
      */
    Topology CurrentParm_;
    /** \brief List of charges that have been converted to Amber units */
    std::vector<double> atom_charge_;
    /** \brief Is our site shape a sphere? If no, it's a box. */
    bool sphere_;
    /** \brief Data set list for all data sets created here */
    DataSetList myDSL_;
    /** \brief List of each peak location */
    std::vector<Vec3> peaks_;
    /** \brief List of the first atom and last atoms of each solvent residue */
    std::vector<Residue> solvent_residues_;
    /** \brief Total number of frames */
    int Nframes_;
    /** \brief Keep track if our cutoff overflowed our box coordinates... */
    bool overflow_;

    // ------------------- Functions -------------------
    int SetupParms(Topology*);
    double Calculate_Energy(Frame*, Residue const&);

    Action::RetType Init(ArgList&, TopologyList*, DataSetList*, DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);
    void Print();
    // Custom Do- routines
    Action::RetType DoPureWater(int, Frame*);
    Action::RetType DoSPAM(int, Frame*);
};

inline bool inside_box(Vec3 gp, Vec3 pt, double edge);
inline bool inside_sphere(Vec3 gp, Vec3 pt, double rad2);

#endif
