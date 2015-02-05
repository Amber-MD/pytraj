#include "Action_ReplicateCell.h"
#include "CpptrajStdio.h"
#include "ParmFile.h"

// CONSTRUCTOR
Action_ReplicateCell::Action_ReplicateCell() : 
  coords_(0), ncopies_(0), ensembleNum_(-1) {} 

void Action_ReplicateCell::Help() {
  mprintf("\t[out <traj filename>] [parmout <parm filename>] [name <dsname>]\n"
          "\t{ all | dir <XYZ> [dir <XYZ> ...] } [<mask>]\n"
          "  Replicate unit cell in specified (or all) directions for atoms in <mask>.\n"
          "    <XYZ>: X= 1, 0, -1, replicate in specified direction (e.g. 100 is +X only)\n");
}

// Action_ReplicateCell::Init()
Action::RetType Action_ReplicateCell::Init(ArgList& actionArgs, TopologyList* PFL, DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  // Require imaging.
  image_.InitImaging( true );
  // Set up output traj
  trajfilename_ = actionArgs.GetStringKey("out");
  parmfilename_ = actionArgs.GetStringKey("parmout");
  Topology* tempParm = PFL->GetParm( actionArgs );
  bool setAll = actionArgs.hasKey("all");
  std::string dsname = actionArgs.GetStringKey("name");
  if (!dsname.empty()) {
    coords_ = (DataSet_Coords*)DSL->AddSet(DataSet::COORDS, dsname, "RCELL");
    if (coords_ == 0) return Action::ERR;
  }
  if (trajfilename_.empty() && coords_ == 0) {
    mprinterr("Error: Either 'out <traj filename> or 'name <dsname>' must be specified.\n");
    return Action::ERR;
  }
  // Get Mask
  Mask1_.SetMaskString( actionArgs.GetMaskNext() );

  // Determine which directions to set
  if (setAll) {
    for (int ix = -1; ix < 2; ix++)
      for (int iy = -1; iy < 2; iy++)
        for (int iz = -1; iz < 2; iz++) {
          directionArray_.push_back( ix );
          directionArray_.push_back( iy );
          directionArray_.push_back( iz );
        }
  } else {
    std::string dirstring = actionArgs.GetStringKey("dir");
    while (!dirstring.empty()) {
      std::vector<int> ixyz(3, -2);
      std::vector<int>::iterator iptr = ixyz.begin();
      for (std::string::const_iterator c = dirstring.begin();
                                       c != dirstring.end(); ++c)
      {
        if (iptr == ixyz.end()) {
          mprinterr("Error: 'dir' string has too many characters.\n");
          return Action::ERR;
        }
        int sign = 1;
        if      (*c == '+') ++c;
        else if (*c == '-') { sign = -1; ++c; }
        
        if      (*c == '1') *iptr = 1 * sign;
        else if (*c == '0') *iptr = 0;
        ++iptr;
      }
      mprintf("DEBUG: %s = %i %i %i\n", dirstring.c_str(), ixyz[0], ixyz[1], ixyz[2]);
      directionArray_.push_back( ixyz[0] );
      directionArray_.push_back( ixyz[1] );
      directionArray_.push_back( ixyz[2] );
      dirstring = actionArgs.GetStringKey("dir");
    }
  }
  ncopies_ = (int)(directionArray_.size() / 3);
  if (ncopies_ < 1) {
    mprinterr("Error: No directions (or 'all') specified.\n");
    return Action::ERR;
  }
  // Set up output trajectory
  if (!trajfilename_.empty()) {
    if (tempParm == 0) {
      mprinterr("Error: Could not get topology for %s\n", trajfilename_.c_str());
      return Action::ERR;
    }
    outtraj_.SetDebug( debugIn );
    // Initialize output trajectory with remaining arguments
    trajArgs_ = actionArgs.RemainingArgs();
    ensembleNum_ = DSL->EnsembleNum();
  }

  mprintf("    REPLICATE CELL: Replicating cell in %i directions:\n", ncopies_);
  mprintf("\t\t X  Y  Z\n");
  for (unsigned int i = 0; i != directionArray_.size(); i += 3)
    mprintf("\t\t%2i %2i %2i\n", directionArray_[i], 
            directionArray_[i+1], directionArray_[i+2]);
  mprintf("\tUsing atoms in mask '%s'\n", Mask1_.MaskString());
  if (!trajfilename_.empty())
    mprintf("\tWriting to trajectory %s\n", trajfilename_.c_str());
  if (!parmfilename_.empty())
    mprintf("\tWriting topology %s\n", parmfilename_.c_str());
  if (coords_ != 0)
    mprintf("\tSaving coords to data set %s\n", coords_->Legend().c_str());

  return Action::OK;
}

// Action_ReplicateCell::Setup()
/** Determine what atoms each mask pertains to for the current parm file.
  */
Action::RetType Action_ReplicateCell::Setup(Topology* currentParm, Topology** parmAddress) {
  if (currentParm->SetupIntegerMask( Mask1_ )) return Action::ERR;
  mprintf("\t%s (%i atoms)\n",Mask1_.MaskString(), Mask1_.Nselected());
  if (Mask1_.None()) {
    mprintf("Warning: One or both masks have no atoms.\n");
    return Action::ERR;
  }
  // Set up imaging info for this parm
  image_.SetupImaging( currentParm->BoxType() );
  if (!image_.ImagingEnabled()) {
    mprintf("Warning: Imaging cannot be performed for topology %s\n", currentParm->c_str());
    return Action::ERR;
  }
  // Create combined topology.
  if (combinedTop_.Natom() > 0) {
    // Topology already set up. Check that # atoms matches.
    if (Mask1_.Nselected() * ncopies_ != combinedTop_.Natom()) {
      mprintf("Warning: Unit cell can currently only be replicated for"
              " topologies with same # atoms.\n");
      return Action::ERR;
    }
    // Otherwise assume top does not change.
  } else {
    // Set up topology and frame.
    Topology* stripParm = currentParm->modifyStateByMask( Mask1_ );
    if (stripParm == 0) return Action::ERR;
    for (int cell = 0; cell != ncopies_; cell++)
      combinedTop_.AppendTop( *stripParm );
    combinedTop_.Brief("Combined parm:");
    delete stripParm;
    if (!parmfilename_.empty()) {
      ParmFile pfile;
      if (pfile.WriteTopology(combinedTop_, parmfilename_, ParmFile::UNKNOWN_PARM, 0)) {
        mprinterr("Error: Topology file %s not written.\n", parmfilename_.c_str());
        return Action::ERR;
      }
    }
    combinedFrame_.SetupFrameV(combinedTop_.Atoms(), combinedTop_.HasVelInfo(),
                               combinedTop_.NrepDim());
    // Set up COORDS / output traj if necessary.
    if (coords_ != 0)
      coords_->SetTopology( combinedTop_ );
    if (!trajfilename_.empty()) {
      if ( outtraj_.InitEnsembleTrajWrite(trajfilename_, trajArgs_,
                                          &combinedTop_, TrajectoryFile::UNKNOWN_TRAJ,
                                          ensembleNum_) )
        return Action::ERR;
    }
  }

  return Action::OK;
}

// Action_ReplicateCell::DoAction()
Action::RetType Action_ReplicateCell::DoAction(int frameNum, Frame* currentFrame,
                                               Frame** frameAddress)
{
  int idx, newFrameIdx;
  unsigned int id;
  Vec3 frac, t2;

  currentFrame->BoxCrd().ToRecip(ucell_, recip_);
  int shift = Mask1_.Nselected() * 3;
# ifdef _OPENMP
# pragma omp parallel private(idx, newFrameIdx, id) firstprivate(frac, t2)
  {
# pragma omp for
# endif
  for (idx = 0; idx < Mask1_.Nselected(); idx++) {
    // Convert to fractional coords
    frac = recip_ * Vec3(currentFrame->XYZ( Mask1_[idx] ));
    // replicate in each direction
    newFrameIdx = idx * 3;
    for (id = 0; id != directionArray_.size(); id+=3, newFrameIdx += shift)
    {
       // Convert back to Cartesian coords.
       t2 = ucell_.TransposeMult(frac + Vec3(directionArray_[id  ],
                                             directionArray_[id+1],
                                             directionArray_[id+2]));
       combinedFrame_[newFrameIdx  ] = t2[0];
       combinedFrame_[newFrameIdx+1] = t2[1];
       combinedFrame_[newFrameIdx+2] = t2[2];
    }
  }
# ifdef _OPENMP
  }
# endif
  if (!trajfilename_.empty()) {
    if (outtraj_.WriteFrame(frameNum, &combinedTop_, combinedFrame_)!=0)
      return Action::ERR;
  }
  if (coords_ != 0)
    coords_->AddFrame( combinedFrame_ );
  
  return Action::OK;
}
