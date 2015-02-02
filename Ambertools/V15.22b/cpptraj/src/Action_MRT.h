#ifndef INC_ACTION_MRT_H
#define INC_ACTION_MRT_H
#include "Action.h"
/** Perform mean residence time calculations.
  * \author Original code by: Hannes H. Loeffler
  *         2005-2007: Lab. of Molecular Design, ICMB, the University of Tokyo, Japan
  *         2008-2010: STFC Daresbury Laboratory, UK
  * \author Adapted to C++ by Daniel R. Roe.
  */
class Action_MRT : public Action {
  public:
    Action_MRT();

    static DispatchObject* Alloc() { return (DispatchObject*)new Action_MRT(); }
    static void Help();

  private:
    Action::RetType Init(ArgList&, TopologyList*, FrameList*, DataSetList*,
                          DataFileList*, int);

    CpptrajFile outfile_;
    double time_; // darg1
    int nStar_;              ///< Nsteps not counted as having left/entered
    double lowerCutoff2_;    ///< define lower cutoff for inside region
    double upperCutoff2_;    ///< define upper cutoff for inside region
    std::string autoCorr_;   ///< filename for acf output
    int wSize_;              ///< window size in steps
    int nOffset_;            ///< offset between windows
    int idxMaxWin_;          ///< maximum number of parallel windows
    AtomMask solventmask_;   ///< Solvent mask expression
    AtomMask solutemask_;    ///< Solute mask expression
};
#endif
