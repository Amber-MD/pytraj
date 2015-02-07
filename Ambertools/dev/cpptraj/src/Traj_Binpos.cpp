#include "Traj_Binpos.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Traj_Binpos::Traj_Binpos() :
  bpnatom_(0),
  bpnatom3_(0),
  frameSize_(0),
  bpbuffer_(0)
{}

Traj_Binpos::~Traj_Binpos() {
  if (bpbuffer_!=0) delete[] bpbuffer_;
}

bool Traj_Binpos::ID_TrajFormat(CpptrajFile& fileIn) {
  unsigned char buffer[4];
  buffer[0] = ' ';
  buffer[1] = ' ';
  buffer[2] = ' ';
  buffer[3] = ' ';
  if ( fileIn.OpenFile() ) return false;
  fileIn.Read(buffer, 4);
  fileIn.CloseFile();
  // Check for the magic header of the Scripps binary format.
  if (buffer[0] == 'f' &&
      buffer[1] == 'x' &&
      buffer[2] == 'y' &&
      buffer[3] == 'z')
    return true;
  return false;
}

int Traj_Binpos::openTrajin() {
  unsigned char buffer[4]; 
  if (file_.OpenFile()) return 1;
  // Read past magic header
  if (file_.Read(buffer, 4) != 4) return 1;
  return 0;
}

void Traj_Binpos::closeTraj() {
  file_.CloseFile();
}

int Traj_Binpos::setupTrajin(std::string const& fname, Topology* trajParm)
{
  if (file_.SetupRead( fname, 0 )) return TRAJIN_ERR;
  // Open - reads past the 4 byte header
  if (openTrajin()) return TRAJIN_ERR;
  // Binpos file is 4 byte header followed by binpos frame records.
  // Each frame record is an int (# atoms) followed by 3*natom
  // floats (the xyz coords).
  // First assume binpos file has consistent # of atoms in each record
  // and try to determine # of frames.
  file_.Read(&bpnatom_, sizeof(int));
  if (bpnatom_ != trajParm->Natom()) {
    mprinterr("Error: # of atoms in binpos file frame 1 (%i) is not equal to\n",bpnatom_);
    mprinterr("Error: the # of atoms in associated parm %s (%i)\n",
              trajParm->c_str(), trajParm->Natom());
    return TRAJIN_ERR;
  }
  bpnatom3_ = bpnatom_ * 3;
  frameSize_ = (size_t)bpnatom3_ * sizeof(float);
  off_t framesize = (off_t)frameSize_ + sizeof(int);
  off_t filesize = file_.UncompressedSize();
  int Frames = 0;
  if (filesize < 1) {
    mprintf("Warning: binpos: Could not determine file size for # frames prediction.\n");
    mprintf("Warning: This is normal for bzip2 files.\n");
    Frames = TRAJIN_UNK;
  } else { 
    filesize -= 4; // Subtract 4 byte header
    Frames = (int)(filesize / framesize);
    if ( (filesize % framesize) != 0 ) {
      mprintf("Warning: %s: Could not accurately predict # frames. This usually\n",
              "Warning:  indicates a corrupted trajectory or topology/trajectory\n"
              "Warning:  mismatch. Will attempt to read %i frames.\n",
              file_.Filename().base(), Frames);
    }
  }
  mprintf("\t%i atoms, framesize=%lu, filesize=%lu, #Frames=%i\n", 
          bpnatom_, framesize, filesize, Frames);
  // Allocate space to read in floats
  if (bpbuffer_!=0)
    delete[] bpbuffer_;
  bpbuffer_ = new float[ bpnatom3_ ];
  closeTraj();
  return Frames;
}

int Traj_Binpos::readFrame(int set, Frame& frameIn) {
  int natoms;
  // Seek
  file_.Seek( ((off_t)set * (frameSize_ + sizeof(int))) + 4 );
  // Read past natom
  if (file_.Read(&natoms, sizeof(int))<1)
    return 1;
  // Sanity check
  if ( natoms != bpnatom_ ) {
    mprinterr("Error: Reading of binpos files with varying # of atoms is not supported.\n");
    return 1;
  }
  // Read coords
  file_.Read(bpbuffer_, frameSize_);
  // Convert float to double
  for (int i = 0; i < bpnatom3_; ++i)
    frameIn[i] = (double)bpbuffer_[i];
  return 0;
}

int Traj_Binpos::setupTrajout(std::string const& fname, Topology* trajParm,
                              CoordinateInfo const& cInfoIn,
                              int NframesToWrite, bool append)
{
  if (!append) {
    if (file_.SetupWrite( fname, debug_ )) return 1;
    unsigned char buffer[4];
    bpnatom_ = trajParm->Natom();
    bpnatom3_ = bpnatom_ * 3;
    frameSize_ = (size_t)bpnatom3_ * sizeof(float);
    // Allocate space to write out floats
    if (bpbuffer_!=0)
      delete[] bpbuffer_;
    bpbuffer_ = new float[ bpnatom3_ ];
    if (CoordInfo().HasBox()) 
      mprintf("Warning: BINPOS format does not support writing of box coordinates.\n");
    if (file_.OpenFile()) return 1;
    // Always write header
    buffer[0] = 'f';
    buffer[1] = 'x';
    buffer[2] = 'y';
    buffer[3] = 'z';
    file_.Write(buffer, 4);
  } else {
    if ( setupTrajin( fname, trajParm ) == TRAJIN_ERR ) return 1;
    // Re-open for appending
    if (file_.SetupAppend( fname, debug_ )) return 1;
    if (file_.OpenFile()) return 1;
  }
  return 0;
}

int Traj_Binpos::writeFrame(int set, Frame const& frameOut) {
  file_.Write( &bpnatom_, sizeof(int) );
  // Convert double to float
  for (int i = 0; i < bpnatom3_; ++i)
    bpbuffer_[i] = (float)frameOut[i];
  if (file_.Write( bpbuffer_, frameSize_ )) return 1;
  return 0;
}

void Traj_Binpos::Info() {
  mprintf("is a BINPOS file");
}
 
