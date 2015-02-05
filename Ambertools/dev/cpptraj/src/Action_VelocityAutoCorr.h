#ifndef INC_ACTION_VELOCITYAUTOCORR_H
#define INC_ACTION_VELOCITYAUTOCORR_H
#include "Action.h"
#include "DataSet_Vector.h"
class Action_VelocityAutoCorr : public Action {
  public:
    Action_VelocityAutoCorr();
    static DispatchObject* Alloc() { return (DispatchObject*)new Action_VelocityAutoCorr(); }
    static void Help();
  private:
    Action::RetType Init(ArgList&, TopologyList*, DataSetList*, DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);
    void Print();

    bool useVelInfo_;     ///< If true use actual velocities in frame if present
    bool useFFT_;         ///< Use FFT to calculate VAC functions
    bool normalize_;      ///< Normalize VAC fn to 1.0 
    AtomMask mask_;       ///< Atoms to calculate VAC fn for.
    Frame previousFrame_; ///< Hold previous frame coords (!useVelInfo only)
    typedef DataSet_Vector Varray;
    typedef std::vector<Varray> VelArray;
    VelArray Vel_; ///< Hold velocity info for each selected atom at each frame.
    DataSet* VAC_; ///< Hold values of the velocity auto-correlation function
    double tstep_; ///< Time step
    int maxLag_;   ///< Maximum lag to calculate VAC fn for.
};
#endif
