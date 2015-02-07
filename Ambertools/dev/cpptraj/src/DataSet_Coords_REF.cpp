#include "DataSet_Coords_REF.h"
#include "CpptrajStdio.h"
#include "Trajin_Single.h"

void DataSet_Coords_REF::Info() const {
  //if (!tag_.empty())
  //  mprintf(" %s", tag_.c_str());
  if (!name_.empty() && Name() != name_.Full())
    mprintf(" '%s'", name_.full());
  if (num_ > -1) mprintf(", refindex %i", num_);
  mprintf(" (%i atoms)", Top().Natom());
}

int DataSet_Coords_REF::LoadRef(std::string const& fname, Topology const& parmIn, int dbg)
{
  ArgList blank;
  return SetupRefFrame(fname, "", parmIn, blank, -1);
}

int DataSet_Coords_REF::SetupRefFrame(std::string const& fname, std::string const& nameIn,
                                      Topology const& parmIn, ArgList& argIn, int refidx)
{
  // Set up trajectory - false = do not modify box info
  Trajin_Single traj;
  //traj.SetDebug( debug_ );
  if ( traj.SetupTrajRead( fname, argIn, (Topology*)&parmIn, false ) ) { // FIXME: Fix cast
    mprinterr("Error: reference: Could not set up trajectory.\n");
    return 1;
  }
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
  // Set up reference frame
  if (frame_.SetupFrameV(parmIn.Atoms(), traj.TrajCoordInfo()))
    return 1;
  // Read reference frame
  traj.ReadTrajFrame( traj.Start(), frame_ );
  traj.EndTraj();
  name_ = traj.TrajFilename();
  num_ = refidx;
  SetTopology( parmIn );
  // If default name is empty use full trajectory file name.
  std::string setname;
  if (nameIn.empty())
    setname = traj.TrajFilename().Full();
  else
    setname = nameIn;
  SetupSet( setname, traj.Start()+1, "" );
  SetLegend( traj.Title() );
  return 0;
}

/** Currently used by 'reference' and 'atommap' */
int DataSet_Coords_REF::StripRef(AtomMask const& stripMask) {
    Frame stripFrame( frame_, stripMask );
    Topology* stripParm = Top().modifyStateByMask( stripMask );
    if (stripParm == 0) {
      mprinterr("Error: Could not create stripped reference topology.\n");
      return 1;
    }
    stripParm->Brief("Stripped ref parm:");
    frame_ = stripFrame;
    SetTopology( *stripParm );
    delete stripParm; // OK to free, parm has been copied by SetTopology.
    return 0;
}
