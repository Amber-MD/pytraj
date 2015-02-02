// TrajoutList
#include "TrajoutList.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // integerToString 
//#incl ude "MpiRoutines.h" //worldsize

TrajoutList::TrajoutList() : debug_(0) { }

TrajoutList::~TrajoutList() {
  Clear();
}

void TrajoutList::SetDebug(int debugIn) {
  debug_ = debugIn;
  if (debug_ > 0)
    mprintf("TrajoutList debug level set to %i\n", debug_);
}

void TrajoutList::Clear() {
  for (ListType::iterator traj = trajout_.begin(); traj != trajout_.end(); ++traj) 
    delete *traj;
  trajout_.clear();
  trajoutArgs_.clear();
}

// TrajoutList::AddEnsembleTrajout()
int TrajoutList::AddEnsembleTrajout(ArgList const& argIn, TopologyList const& topListIn,
                                    int member)
{
  // Make a copy of input arg list so that args remain unmarked for the next
  // member of the ensemble.
  ArgList args = argIn;
  std::string filename = args.GetStringNext();
  if (filename.empty()) {
    mprinterr("Error: TrajoutList::AddEnsemble: Called with null filename.\n");
    return 1;
  }
  std::string onlyMembers = args.GetStringKey("onlymembers");
  if (!onlyMembers.empty()) {
    // Range of ensemble members to write for. If this member is
    // not on the list, exit.
    Range members( onlyMembers );
    if (members.Empty()) {
      mprinterr("Error: onlymembers: Invalid range (%s)\n", onlyMembers.c_str());
      return 1;
    }
    bool is_a_member = false;
    for (Range::const_iterator mstr = members.begin(); mstr != members.end(); ++mstr)
      if ( member == *mstr ) {
        is_a_member = true;
        break;
      }
    if ( !is_a_member ) {
      mprintf("trajout %s: Member %i is not on onlymembers list '%s'; will not be written to.\n",
              filename.c_str(), member, onlyMembers.c_str());
      return 0;
    }
  }
  // Trajectory format keywords take precendence over extension.
  TrajectoryFile::TrajFormatType outtrajFmt = TrajectoryFile::GetFormatFromArg(args);
  if (outtrajFmt == TrajectoryFile::AMBERTRAJ) {
    // Before filename is modified see if extension is recognized
    FileName tempName;
    tempName.SetFileName( filename );
    outtrajFmt = TrajectoryFile::GetTypeFromExtension( tempName.Ext() );
  }
  // Modify filename by member
  filename += ("." + integerToString( member ));
  return AddTrajout( filename, args, topListIn, outtrajFmt );
}

// TrajoutList::AddTrajout()
int TrajoutList::AddTrajout(ArgList const& argIn, TopologyList const& topListIn) {
  // Since we need to check if this filename is in use in order to prevent
  // overwrites, determine the filename here.
  ArgList args = argIn;
  std::string filename = args.GetStringNext();
  if (filename.empty()) {
    mprinterr("Internal Error: TrajoutList::Add() called with empty filename.\n");
    return 1;
  }
  int err = AddTrajout( filename, args, topListIn, TrajectoryFile::UNKNOWN_TRAJ );
  // For setting up ensemble later, save trajout arg.
  if (err == 0) trajoutArgs_.push_back( argIn );
  return err;
}

// TrajoutList::AddTrajout()
/** Add trajectory to the trajectory list as an output trajectory. 
  * Associate the trajectory with one of the parm files in the 
  * TopologyList. 
  */
int TrajoutList::AddTrajout(std::string const& filename, ArgList& argIn, 
                            TopologyList const& topListIn,
                            TrajectoryFile::TrajFormatType fmtIn) 
{
  // Check if filename is in use
  for (ListType::const_iterator to = trajout_.begin();
                                to != trajout_.end(); ++to)
  {
    if ( (*to)->TrajFilename().Full() == filename ) { 
      mprinterr("Error: trajout: Filename %s already in use.\n",filename.c_str());
      return 1;
    }
  }
  // Create trajout.
  Trajout *traj = new Trajout();
  if (traj==0) {
    mprinterr("Error: TrajoutList::Add: Could not allocate memory for traj.\n");
    return 1;
  }
  // Get parm from TopologyList based on args
  Topology* tempParm = topListIn.GetParm( argIn );
  traj->SetDebug(debug_);
  // Default to AMBERTRAJ; format can be changed via args in the arg list
  if (traj->InitTrajWrite(filename, argIn, tempParm, fmtIn)) {
    mprinterr("Error: trajout: Could not set up trajectory.\n");
    delete traj;
    return 1;
  }

  // Add to trajectory file list
  trajout_.push_back(traj);

  return 0;
}

// TrajoutList::WriteTrajout()
/** Go through each output traj, call write. The first time CurrentParm
  * matches the parm the trajectory was originally set up with it will
  * be opened, no need to call BeginTraj.
  */ 
int TrajoutList::WriteTrajout(int set, Topology *CurrentParm, Frame *CurrentFrame) { 
  for (ListType::iterator traj = trajout_.begin(); traj != trajout_.end(); ++traj) 
  {
    if ( (*traj)->WriteFrame(set, CurrentParm, *CurrentFrame) ) {
      mprinterr("Error writing output trajectory, frame %i.\n", set+1);
      return 1;
    }
  }
  return 0;
}

// TrajoutList::CloseTrajout()
/** Close output trajectories. Called after input traj processing completed. */
void TrajoutList::CloseTrajout() {
  for (ListType::iterator traj = trajout_.begin(); traj != trajout_.end(); ++traj)
    (*traj)->EndTraj();
  Clear();
}

// TrajoutList::List()
void TrajoutList::List() const {
  if (!trajout_.empty()) {
    mprintf("\nOUTPUT TRAJECTORIES:\n");
    for (ListType::const_iterator traj = trajout_.begin(); traj != trajout_.end(); ++traj)
      (*traj)->PrintInfo( 1 );
  }
}
