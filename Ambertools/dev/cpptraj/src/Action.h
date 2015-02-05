#ifndef INC_ACTION_H
#define INC_ACTION_H
#include "DispatchObject.h"
#include "ArgList.h"
#include "DataFileList.h"
#include "DataSetList.h"
#include "TopologyList.h"
// Class: Action 
/// The abstract base class that all other actions inherit. 
/** By convention actions have 3 main phases: Init, Setup, and DoAction.
  * Init is used to initialize the action, make sure that all arguments
  * for the action are correct, and add any DataSets/DataFiles which will
  * be used by the action. Setup will set up the action for a specific
  * Topology file. DoAction will perform the action on a given frame.
  * A fourth function, Print, is for any additional calculations or output 
  * the action may require once all frames are processed.
  * Note that Setup and DoAction take Topology and Frame pointer addresses 
  * respectively as their final arguments; this allows actions to modify
  * the current Topology/Frame by replacing them with ones inside the Action
  * itself. This is more memory-hungry but allows the Topology/Frame to be
  * manipulated quickly, and is also easier to return to the original
  * Topology/Frame via USEORIGINALFRAME (see e.g. Action_Unstrip in 
  * Action_Strip.cpp).
  */
class Action : public DispatchObject {
  public:
    /// Enumerate potential return states from Init, Setup, and DoAction.
    enum RetType { OK=0, ///< Everything OK, normal return.
                   ERR,  ///< Problem occurred.
                   USEORIGINALFRAME, ///< Return to unmodified frame/topology.
                   SUPPRESSCOORDOUTPUT ///< Skip remaining actions and traj output.
    };
    /// Destructor - virtual since this class is inherited
    virtual ~Action() {}
    /// Initialize action
    /** Process input args, set up any DataSets or DataFiles, set debug level */
    virtual RetType Init(ArgList&, TopologyList*, DataSetList*, DataFileList*, int) = 0;
    /// Set up action for given Topology
    virtual RetType Setup(Topology*,Topology**) = 0;
    /// Perform action for given frame number and Frame.
    virtual RetType DoAction(int,Frame*,Frame**) = 0;
    /// Print anything besides datasets, called at end of execution
    /** Perform any output not related to master dataset output, or any 
      * necessary post-trajectory processing calculations.
      */
    virtual void Print() = 0;
};
#endif  
