#include "Traj_Tinker.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Traj_Tinker::Traj_Tinker() : tinkerTop_(0), currentSet_(0) {}

bool Traj_Tinker::ID_TrajFormat(CpptrajFile& fileIn) {
  return TinkerFile::ID_Tinker( fileIn );
}

// Traj_Tinker::openTrajin()
int Traj_Tinker::openTrajin() {
  currentSet_ = 0;
  return file_.OpenTinker();
}

// Traj_Tinker::closeTraj() {
void Traj_Tinker::closeTraj() {
  file_.CloseFile();
}

// Traj_Tinker::setupTrajin()
/** See how many XYZ records are in file, make sure num atoms match
  * parm and each frame.
  */
int Traj_Tinker::setupTrajin(std::string const& fname, Topology* trajParm)
{
  int Frames=0;
  file_.SetTinkerName( fname );
  if (file_.OpenTinker()) return TRAJIN_ERR;
  // Check # atoms in tinker file against # atoms in parm
  if (file_.TinkerNatom() != trajParm->Natom()) {
    mprinterr("Error: Number of atoms in Tinker file %s (%i) does not\n",
              file_.Filename().base(), file_.TinkerNatom());
    mprinterr("Error: match number in associated parmtop (%i)!\n",trajParm->Natom());
    return TRAJIN_ERR;
  }
  // Set title
  SetTitle( file_.TinkerTitle() );
  // Set traj info - no velocity, temperature, time.
  SetCoordInfo( CoordinateInfo(file_.TinkerBox(), false, false, false) );
  // Count how many frames can be read.
  int readMoreFrames = 1;
  while (readMoreFrames == 1) {
    readMoreFrames = file_.NextTinkerFrame();
    Frames += readMoreFrames;
  }
  if (readMoreFrames == -1) {
    mprintf("Warning: An error occurred while determining number of frames in Tinker file.\n"
            "Warning: Will attempt to read %i frames.\n", Frames);
  }
      
  file_.CloseFile();
  //if (debug_>0)
    mprintf("\tTinker file %s has %i frames.\n",file_.Filename().base(), Frames);

  return Frames;
}

// Traj_Tinker::readFrame()
int Traj_Tinker::readFrame(int set, Frame& frameIn) {
  if (set < currentSet_) {
    file_.Rewind();
    currentSet_ = 0;
  }
  // Position file at specified set
  while (currentSet_ < set) {
    if (file_.NextTinkerFrame()==-1) return 1;
    currentSet_++;
  }
  if (file_.ReadNextTinkerFrame( frameIn.xAddress(), frameIn.bAddress() ) != 1)
    return 1;
  currentSet_++;
  return 0;
}

// Traj_Tinker::setupTrajout()
/** Set parm information required for write, and check write mode against
  * number of frames to be written.
  */
int Traj_Tinker::setupTrajout(std::string const& fname, Topology* trajParm,
                              CoordinateInfo const& cInfoIn,
                              int NframesToWrite, bool append)
{
  // No write yet.
  mprinterr("Error: WRITE not yet implemented for Tinker files.\n");
  return 1;
/*  if (trajParm==0) return 1;
  mol2Top_ = trajParm;
  // Set up file
  if (append && mol2WriteMode_ != MULTI) {
    if (file_.SetupAppend( fname, debug_)) return 1;
  } else {
    if (mol2WriteMode_ == MULTI)
      mprintf("Warning: 'append' not compatible with 'multi' mol2 write.\n");
    if (file_.SetupWrite( fname, debug_ )) return 1;
  }
  // If writing more than 1 frame and not writing 1 pdb per frame, 
  // use @<TRIPOS>MOLECULE keyword to separate frames.
  if (append || (mol2WriteMode_==SINGLE && NframesToWrite>1)) 
    mol2WriteMode_ = MOL;
  // Set # atoms; if more than 99999 atoms the file may not write correctly
  file_.SetMol2Natoms( mol2Top_->Natom() );
  if (file_.Mol2Natoms() > 99999) {
    mprintf("Warning: %s: Large # of atoms (%i > 99999) for mol2 format.\n",
            file_.Filename().base(), file_.Mol2Natoms());
    mprintf("Warning: File may not write correctly.\n");
  }
  // TODO: Change this, right now for backwards compat. only!
  // If all charges == 0 set noCharges.
  hasCharges_ = false;
  for (Topology::atom_iterator atom = mol2Top_->begin(); atom != mol2Top_->end(); atom++)
  {
    if ( (*atom).Charge() != 0 ) {
      hasCharges_ = true;
      break;
    }
  }
  // Set Title
  if (Title().empty())
    SetTitle("Cpptraj generated mol2 file.");
  file_.SetMol2Title( Title() );
  // Set up number of bonds
  file_.SetMol2Nbonds( mol2Top_->Bonds().size() + mol2Top_->BondsH().size() );
  // Open here if writing to a single file
  if ( mol2WriteMode_ != MULTI )
    return file_.OpenFile(); 
  return 0;*/
}

// Traj_Tinker::writeFrame()
int Traj_Tinker::writeFrame(int set, Frame const& frameOut) {
  return 1;
}
 
// Traj_Tinker::info()
void Traj_Tinker::Info() {
  mprintf("is a Tinker file");
}
