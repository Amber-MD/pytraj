#ifndef INC_ACTION_MASK_H
#define INC_ACTION_MASK_H
#include "Action.h"
#include "TrajectoryFile.h"
#include "DataSet.h"
// Class: Action_Mask_2
/// Print out all atoms selected by a mask for each frame.
/** This allows use of distance-dependent masks. This does NOT modify the
  * frame or parm. 
  */
class Action_Mask_2: public Action {
  public:
    Action_Mask_2();
    static DispatchObject* Alloc() { return (DispatchObject*)new Action_Mask_2(); }
    static void Help();
  private:
    Action::RetType Init(ArgList&, TopologyList*, DataSetList*, DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);
    void Print();

    AtomMask Mask1_;         ///< Atoms which will be selected each frame
    CpptrajFile outfile_;    ///< File to write selected atom info to
    std::string maskpdb_;    ///< Traj output file name
    Topology* CurrentParm_;
    int debug_;
    TrajectoryFile::TrajFormatType trajFmt_; ///< Output trajectory format
    const char* trajOpt_;    ///< Output trajectory options

    DataSet* atomDs_;        /// holding atom index for pytraj
    DataSet* frameIndexDs_;  /// holding frame index for pytraj
};
#endif  
