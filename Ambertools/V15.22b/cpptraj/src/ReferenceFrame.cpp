#include "ReferenceFrame.h"
#include "Trajin_Single.h"
#include "CpptrajStdio.h"

// DESTRUCTOR
// Memory must be freed with ClearRef();
ReferenceFrame::~ReferenceFrame() {}

void ReferenceFrame::ClearRef() {
  if (frame_ != 0) delete frame_;
  if (parm_ != 0) delete parm_;
  frame_ = 0;
  parm_ = 0;
  num_ = 0;
  tag_.clear();
  name_.clear();
}

int ReferenceFrame::LoadRef(std::string const& fname, Topology* parmIn, int debugIn)
{
  ArgList argIn;
  return LoadRef( fname, argIn, parmIn, std::string(), debugIn);
}

// ReferenceFrame::LoadRef()
int ReferenceFrame::LoadRef(std::string const& fname, ArgList& argIn, 
                            Topology* parmIn, std::string const& maskexpr,
                            int debugIn)
{
  Trajin_Single traj;
  traj.SetDebug(debugIn);
  num_ = -1; // This is -1 as long as ref is not fully set up.
  // Set up trajectory - false = do not modify box info
  if ( traj.SetupTrajRead( fname, argIn, parmIn, false ) ) {
    mprinterr("Error: reference: Could not set up trajectory.\n");
    return 1;
  }
  // Check for tag - done after SetupTrajRead so traj can process args
  tag_ = argIn.getNextTag();
  // Check number of frames to be read
  int trajFrames = traj.TotalReadFrames();
  if (trajFrames < 1) {
    mprinterr("Error: No frames could be read for reference '%s'\n", traj.TrajFilename().full());
    return 1;
  } else if (trajFrames > 1)
    mprintf("Warning: Reference has %i frames, only reading frame %i\n", 
            trajFrames, traj.Start()+1);
  // Start trajectory read
  if ( traj.BeginTraj(false) ) {
    mprinterr("Error: Could not open reference '%s'\n.", traj.TrajFilename().full());
    return 1;
  }
  mprintf("\t");
  traj.PrintInfo(1);
  // Set topology; remove old topology if it was set.
  if (parm_ != 0) delete parm_;
  parm_ = new Topology(*parmIn);
  // Set up input frame
  if (frame_ == 0) frame_ = new Frame();
  frame_->SetupFrameV( parm_->Atoms(), traj.HasVelocity(), traj.NreplicaDimension() );
  // Read reference frame
  traj.ReadTrajFrame( traj.Start(), *frame_ );
  // If a mask expression was specified, strip to match the expression.
  if (!maskexpr.empty()) {
    AtomMask stripMask( maskexpr );
    if (parm_->SetupIntegerMask(stripMask)) return 1;
    if ( StripRef( stripMask ) ) {
      delete frame_;
      frame_ = 0;
      return 1;
    }
  }
  // Set name and number.
  name_ = traj.TrajFilename();
  num_ = traj.Start();
  // Set topology reference coords.
  parm_->SetReferenceCoords( *frame_ );
  return 0;
}

/** Strip reference to match mask */
int ReferenceFrame::StripRef(AtomMask const& stripMask) {
  if (stripMask.None()) {
    mprinterr("Error: No atoms kept for reference.\n");
    return 1;
  }
  if (frame_ == 0 || frame_->empty()) {
    mprinterr("Internal Error: Attempting to strip empty reference.\n");
    return 1;
  }
  // Create new stripped frame
  Frame *strippedRefFrame = new Frame( *frame_, stripMask );
  mprintf("\tKept %i atoms in reference.\n", strippedRefFrame->Natom());
  // Create new stripped parm
  Topology *strippedRefParm = parm_->modifyStateByMask( stripMask );
  if (strippedRefParm==0) {
    mprinterr("Error: Could not create stripped reference topology.\n");
    return 1;
  }
  strippedRefParm->Brief("Stripped ref parm:");
  delete frame_;
  frame_ = strippedRefFrame;
  delete parm_;
  parm_ = strippedRefParm;
  return 0;
}

void ReferenceFrame::RefInfo() const {
  if (!tag_.empty())
    mprintf(" %s", tag_.c_str());
  mprintf(" '%s', frame %i\n", name_.full(), num_ + 1);
}
