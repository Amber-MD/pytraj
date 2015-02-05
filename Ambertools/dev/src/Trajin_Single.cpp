#include "Trajin_Single.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // fileExists

// CONSTRUCTOR
Trajin_Single::Trajin_Single() :
  trajio_(0),
  velio_(0),
  trajIsOpen_(false)
{}

// DESTRUCTOR
Trajin_Single::~Trajin_Single() {
  if (trajio_!=0) {
    if (trajIsOpen_) EndTraj();
    delete trajio_;
  }
  if (velio_!=0) delete velio_;
}

int Trajin_Single::SetupTrajRead(std::string const& tnameIn, ArgList& argIn, Topology* tparmIn)
{
  return SetupTrajRead(tnameIn, argIn, tparmIn, true);
}

// Trajin_Single::SetupTrajRead()
int Trajin_Single::SetupTrajRead(std::string const& tnameIn, ArgList& argIn, 
                                 Topology* tparmIn, bool checkBox) 
{
  // Require a filename
  if (tnameIn.empty()) {
    mprinterr("Internal Error: Trajin_Single: No filename given.\n");
    return 1;
  }
  // Check and set associated parm file
  if ( SetTrajParm( tparmIn ) ) return 1;
  // Check that file can be opened. 
  if (!fileExists(tnameIn)) return 1; 
  // Detect file format
  TrajFormatType tformat;
  if ( (trajio_ = DetectFormat( tnameIn, tformat )) == 0 ) {
    mprinterr("Error: Could not determine trajectory %s format.\n", tnameIn.c_str());
    return 1;
  }
  trajio_->SetDebug( debug_ );
  // Set trajectory filename
  SetTrajFileName( tnameIn, true );
  mprintf("\tReading '%s' as %s\n", TrajFilename().full(), TrajectoryFile::FormatString(tformat));
  // Process format-specific read args
  if (trajio_->processReadArgs( argIn )) return 1;
  // Set up the format for reading and get the number of frames.
  if (SetupTrajIO( tnameIn, *trajio_, argIn )) return 1;
  // Check how many frames will actually be read
  if (setupFrameInfo() == 0) return 1;
  // Check traj box info against parm box info
  if (checkBox) { 
    Box parmBox = tparmIn->ParmBox();
    if (CheckBoxInfo(tparmIn->c_str(), parmBox, trajio_->TrajBox())) return 1;
    tparmIn->SetBox( parmBox );
  }
  // Check if a separate mdvel file will be read
  if (argIn.Contains("mdvel")) {
    std::string mdvelname = argIn.GetStringKey("mdvel");
    if (mdvelname.empty()) {
      mprinterr("Error: mdvel: Usage 'mdvel <velocity filename>'\n");
      return 1;
    }
    // Detect mdvel format
    if ( (velio_ = DetectFormat( mdvelname, tformat )) == 0 ) {
      mprinterr("Error: Could not set up velocity file %s for reading.\n",mdvelname.c_str());
      return 1;
    }
    velio_->SetDebug( debug_ );
    // Set up the format for reading mdvel, get # of mdvel frames
    int vel_frames = velio_->setupTrajin(mdvelname, TrajParm());
    if (vel_frames != TotalFrames()) {
      mprinterr("Error: velocity file %s frames (%i) != traj file frames (%i)\n",
                mdvelname.c_str(), vel_frames, TotalFrames());
      return 1;
    }
  }
  return 0;
}

// Trajin_Single::BeginTraj()
int Trajin_Single::BeginTraj(bool showProgress) {
  // Open the trajectory
  if (trajio_->openTrajin()) {
    mprinterr("Error: Trajin_Single::BeginTraj: Could not open %s\n",TrajFilename().base());
    return 1;
  }
  // Open mdvel file if present
  if (velio_!=0 && velio_->openTrajin()) {
    mprinterr("Error: Could not open mdvel file.\n");
    return 1;
  }
  // Set progress bar, start and offset.
  PrepareForRead( showProgress );
  trajIsOpen_ = true;
  return 0;
}

// Trajin_Single::EndTraj()
void Trajin_Single::EndTraj() {
  if (trajIsOpen_) {
    trajio_->closeTraj();
    if (velio_!=0) velio_->closeTraj();
    trajIsOpen_ = false;
  }
}

// Trajin_Single::ReadTrajFrame()
int Trajin_Single::ReadTrajFrame( int currentFrame, Frame& frameIn ) {
  if (trajio_->readFrame(currentFrame, frameIn))
    return 1;
  if (velio_ != 0 && velio_->readVelocity(currentFrame, frameIn))
    return 1;
  //printf("DEBUG:\t%s:  current=%i  target=%i\n",trajName,currentFrame,targetSet);
  return 0;
}

// Trajin_Single::PrintInfo()
void Trajin_Single::PrintInfo(int showExtended) const {
  mprintf("'%s' ",TrajFilename().base());
  trajio_->Info();
  mprintf(", Parm %s",TrajParm()->c_str());
  if (trajio_->HasBox()) mprintf(" (%s box)", trajio_->TrajBox().TypeName());
  if (showExtended==1) PrintFrameInfo(); 
  if (debug_>0)
    mprintf(", %i atoms, Box %i",TrajParm()->Natom(),(int)trajio_->HasBox());
  mprintf("\n");
  if (velio_!=0) {
    mprintf("\tMDVEL: ");
    velio_->Info();
    mprintf("\n");
  }
}

// Trajin_Single::HasVelocity()
bool Trajin_Single::HasVelocity() const {
  if (trajio_!=0) {
    if (velio_ == 0)
      return trajio_->HasV();
    else
      return true;
  }
  return false;
}

int Trajin_Single::NreplicaDimension() const {
  if (trajio_!=0) 
    return trajio_->ReplicaDimensions().Ndims();
  else
    return 0;
}
