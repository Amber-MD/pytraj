#ifndef INC_ACTION_RADIAL_H
#define INC_ACTION_RADIAL_H
#include "Action.h"
#include "ImagedAction.h"
// Class: Action_Radial
/// Calculate the radial distribution (pair correlation) function.
class Action_Radial: public Action, ImagedAction {
  public:
    Action_Radial();
    ~Action_Radial();
    static DispatchObject* Alloc() { return (DispatchObject*)new Action_Radial(); }
    static void Help();
  private:
    Action::RetType Init(ArgList&, TopologyList*, DataSetList*, DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);
    void Print();

    int* RDF_;                ///< Hold bin counts.
    int** rdf_thread_;        ///< Hold bin count on each thread.
    AtomMask Mask1_;          ///< Atoms to calculate RDF for.
    AtomMask Mask2_;          ///< Optional mask to calc RDF to atoms in Mask1.
    AtomMask OuterMask_;      ///< Mask with the most atoms.
    AtomMask InnerMask_;      ///< Mask with the fewest atoms.
    enum RmodeType { NORMAL=0, NO_INTRAMOL, CENTER1, CENTER2 };
    RmodeType rmode_;         ///< Type of calculation to perform.
    Topology* currentParm_;   ///< Current topology, needed for NO_INTERMOL
    int intramol_distances_;  ///< # of intra-molecular distances for NO_INTERMOL.
    bool useVolume_;          ///< If true normalize based on input volume.
    double volume_;           ///< Hold sum of volume for averaging.
    double maximum2_;         ///< Largest distance squared that can be binned.
    double spacing_;          ///< Bin spacing.
    double one_over_spacing_; ///< 1/spacing, used to avoid man division ops.
    int numBins_;             ///< The number of bins.
    int numthreads_;          ///< Number of threads.
    int numFrames_;           ///< Number of frames for which RDF is calcd.
    double density_;          ///< Particle density (molecules/Ang^3).
    DataSet* Dset_;
    DataSet* intrdf_;
    DataSet* rawrdf_;
    int debug_;
};
#endif  
