#include "DataSet_Coords_TRJ.h"
#include "Trajin_Single.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
DataSet_Coords_TRJ::DataSet_Coords_TRJ() :
  DataSet_Coords(TRAJ),
  Traj_(0), 
  currentTrajNum_(-1),
  globalOffset_(0),
  maxFrames_(0),
  deleteTrajectories_(false)
{}

// DESTRUCTOR
DataSet_Coords_TRJ::~DataSet_Coords_TRJ() {
  if (deleteTrajectories_) {
    for (ListType::const_iterator trj = trajinList_.begin(); trj != trajinList_.end(); ++trj)
      delete *trj;
  }
}

// DataSet_Coords_TRJ::SetTrjTopology()
int DataSet_Coords_TRJ::SetTrjTopology( Topology const& parmIn ) {
  if (trajinList_.empty())
    SetTopology( parmIn );
  else {
    if ( parmIn.Natom() != this->Top().Natom() ) {
      mprinterr("Error: For TRAJ data set currently all trajectories must have same number\n"
                "Error:  of atoms: %s (%i) != %s (%i). Recommended course of action is to\n"
                "Error:  create a trajectory where all frames have been stripped to the same\n"
                "Error:  number of atoms first.\n", parmIn.Natom(), this->Top().Natom());
      return 1;
    }
  }
  return 0;
}

// DataSet_Coords_TRJ::UpdateTrjFrames()
int DataSet_Coords_TRJ::UpdateTrjFrames(int Nframes) {
  if (Nframes > 0)
    maxFrames_ += Nframes;
  else {
    mprinterr("Error: Cannot use trajectories with unknown # of frames as data set.\n");
    return 1;
  }
  return 0;
}

/** Add a new Trajin_Single class. */
int DataSet_Coords_TRJ::AddSingleTrajin(std::string const& fname, ArgList& argIn,
                                        Topology* parmIn) // TODO: const&
{
  if (parmIn==0) return 1;
  if (!trajinList_.empty() && !deleteTrajectories_) {
    mprinterr("Internal Error: This DataSet_Coords_TRJ class set up for copies.\n");
    return 1;
  }
  if (SetTrjTopology(*parmIn)) return 1;
  Trajin_Single* trajin = new Trajin_Single();
  if (trajin->SetupTrajRead(fname, argIn, parmIn)) {
    mprinterr("Error: Could not set up trajectory '%s'\n", fname.c_str());
    return 1;
  }
  if (UpdateTrjFrames( trajin->TotalReadFrames() )) return 1;
  trajinList_.push_back( trajin );
  deleteTrajectories_ = true;
  return 0;
} 

/** Add an existing and already set-up Trajin class; this class will not be
  * deleted upon destruction. Intended for use with TrajinList.
  */
int DataSet_Coords_TRJ::AddInputTraj(Trajin* tIn) {
  if (!trajinList_.empty() && deleteTrajectories_) {
    mprinterr("Internal Error: This DataSet_Coords_TRJ class not set up for copies.\n");
    return 1;
  }
  if (tIn == 0) return 1;
  if (SetTrjTopology( *(tIn->TrajParm()) )) return 1;
  if (UpdateTrjFrames( tIn->TotalReadFrames() )) return 1;
  trajinList_.push_back( tIn );
  deleteTrajectories_ = false;
  return 0;
}

void DataSet_Coords_TRJ::GetFrame(int idx, Frame& fIn) {
# ifdef _OPENMP
# pragma omp critical
  {
# endif
  // Determine which trajectory has the desired index
  globalOffset_ = 0;
  int currentMax = 0;
  int desiredTrajNum = 0;
  for (; desiredTrajNum < (int)trajinList_.size(); ++desiredTrajNum) {
    currentMax += trajinList_[desiredTrajNum]->TotalReadFrames();
    if (idx < currentMax) break;
    globalOffset_ += trajinList_[desiredTrajNum]->TotalReadFrames();
  }
  // If desired traj is different than current, open desired traj
  int err = 0;
  if ( desiredTrajNum != currentTrajNum_ ) {
    if (Traj_ != 0) Traj_->EndTraj();
    Trajin* prevTraj = Traj_;
    currentTrajNum_ = desiredTrajNum;
    Traj_ = trajinList_[currentTrajNum_];
    // Set up frame for reading
    bool parmHasChanged;
    if (prevTraj == 0)
      parmHasChanged = true;
    else
      parmHasChanged = ( prevTraj->TrajParm()->Natom() != 
                            Traj_->TrajParm()->Natom()    );
    if ( parmHasChanged || (readFrame_.HasVelocity() != Traj_->TrajCoordInfo().HasVel()) )
      readFrame_.SetupFrameV(Traj_->TrajParm()->Atoms(), Traj_->TrajCoordInfo());
    if (Traj_->BeginTraj(false)) {
      mprinterr("Error: Could not open trajectory %i '%s'\n", desiredTrajNum,
                Traj_->TrajFilename().full());
      err = 1;
    }
  }
  if (err == 0) {
    // Convert desired index into trajectory internal index
    int internalIdx = ((idx - globalOffset_) * Traj_->Offset()) + Traj_->Start();
    // Read desired index
    // TODO: May need to use readFrame here as well...
    if (Traj_->ReadTrajFrame( internalIdx, fIn )) {
      mprinterr("Error: Could not read '%s' frame %i\n", 
                Traj_->TrajFilename().full(), internalIdx + 1);
    }
  }
# ifdef _OPENMP
  }
# endif
}

void DataSet_Coords_TRJ::GetFrame(int idx, Frame& fIn, AtomMask const& mask) {
  GetFrame( idx, readFrame_ );
  fIn.SetFrame( readFrame_, mask );
}

void DataSet_Coords_TRJ::Info() const {
  if (trajinList_.size() == 1)
    mprintf(" (1 trajectory)");
  else
    mprintf(" (%zu trajectories)", trajinList_.size());
}
