// Action_Outtraj 
#include "Action_Outtraj.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Action_Outtraj::Action_Outtraj() : CurrentParm_(0) {} 

void Action_Outtraj::Help() {
  mprintf("\t<filename> [ trajout args ]\n"
          "\t[maxmin <dataset> min <min> max <max>] ...\n"
          "  Like 'trajout', but coordinate output occurs during actions rather than at the end.\n");
}

// Action_Outtraj::Init()
Action::RetType Action_Outtraj::Init(ArgList& actionArgs, TopologyList* PFL, DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  // Set up output traj
  outtraj_.SetDebug(debugIn);
  std::string trajfilename = actionArgs.GetStringNext();
  if (trajfilename.empty()) {
    mprinterr("Error: No filename given.\nError: Usage: ");
    Help();
    return Action::ERR;
  }
  Topology* tempParm = PFL->GetParm(actionArgs);
  if (tempParm==0) {
    mprinterr("Error: OUTTRAJ: Could not get parm for %s\n",trajfilename.c_str());
    return Action::ERR;
  }
  // If maxmin, get the name of the dataset as well as the max and min values.
  double lastmin = 0.0;
  double lastmax = 0.0;
  while ( actionArgs.Contains("maxmin") ) {
    std::string datasetName = actionArgs.GetStringKey("maxmin");
    if (!datasetName.empty()) {
      DataSet* dset = DSL->GetDataSet(datasetName);
      if (dset==0) {
        mprintf("Error: maxmin: Could not get dataset %s\n",datasetName.c_str());
        return Action::ERR;
      } else {
        // Currently only allow int, float, or double datasets
        if (dset->Type() != DataSet::INTEGER &&
            dset->Type() != DataSet::FLOAT &&
            dset->Type() != DataSet::DOUBLE) 
        {
          mprinterr("Error: maxmin: Only int, float, or double dataset (%s) supported.\n",
                  datasetName.c_str());
          return Action::ERR;
        }
        Dsets_.push_back( (DataSet_1D*)dset );
        Max_.push_back( actionArgs.getKeyDouble("max",lastmax) );
        Min_.push_back( actionArgs.getKeyDouble("min",lastmin) );
        lastmax = Max_.back();
        lastmin = Min_.back();
      }
    } else {
      mprinterr("Error: maxmin Usage: maxmin <setname> max <max> min <min>\n");
      return Action::ERR;
    }
  }
  // Initialize output trajectory with remaining arguments
  if ( outtraj_.InitEnsembleTrajWrite(trajfilename, actionArgs.RemainingArgs(), 
                                      tempParm, TrajectoryFile::UNKNOWN_TRAJ,
                                      DSL->EnsembleNum()) ) 
    return Action::ERR;
  mprintf("    OUTTRAJ:");
  outtraj_.PrintInfo(1);
  for (unsigned int ds = 0; ds < Dsets_.size(); ++ds)
    mprintf("\tmaxmin: Printing trajectory frames based on %g <= %s <= %g\n",
            Min_[ds], Dsets_[ds]->Legend().c_str(), Max_[ds]);

  return Action::OK;
} 

// Action_Outtraj::Setup()
Action::RetType Action_Outtraj::Setup(Topology* currentParm, Topology** parmAddress) {
  CurrentParm_ = currentParm; 
  return Action::OK;
}

// Action_Outtraj::DoAction()
/** If a dataset was specified for maxmin, check if this structure
  * satisfies the criteria; if so, write. Otherwise just write.
  */
Action::RetType Action_Outtraj::DoAction(int frameNum, Frame* currentFrame, Frame** frameAddress) {
  // If dataset defined, check if frame is within max/min
  if (!Dsets_.empty()) {
    for (unsigned int ds = 0; ds < Dsets_.size(); ++ds)
    {
      double dVal = Dsets_[ds]->Dval(frameNum);
      //mprintf("DBG: maxmin[%u]: dVal = %f, min = %f, max = %f\n",ds,dVal,Min_[ds],Max_[ds]);
      // If value from dataset not within min/max, exit now.
      if (dVal < Min_[ds] || dVal > Max_[ds]) return Action::OK;
    }
  }
  if ( outtraj_.WriteFrame(frameNum, CurrentParm_, *currentFrame) != 0 ) 
    return Action::ERR;
  return Action::OK;
}

// Action_Outtraj::Print()
/** Close trajectory. Indicate how many frames were actually written.
  */
void Action_Outtraj::Print() {
  mprintf("  OUTTRAJ: [%s] Wrote %i frames.\n",outtraj_.TrajFilename().base(),
          outtraj_.NumFramesProcessed());
  outtraj_.EndTraj();
}
